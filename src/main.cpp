#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>

#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// defines
#define TIME_INTERVAL 0.02
#define TIME_SPAN 2

#define MAX_ACC 0.2
#define CONVERT_M_2_MILES 2.24

#ifdef CONSERVATIVE
#define MAX_SPEED 40
#define N_POINTS 50 // or TIME_SPAN/TIME_INTERVAL/some constant
#define N_POINTS_FOR_SPLINE 3
#define DISTANCE_THRESHOLD 35
#else
#define MAX_SPEED 49
#define N_POINTS 100 // or TIME_SPAN/TIME_INTERVAL/some constant
#define N_POINTS_FOR_SPLINE 5
#define DISTANCE_THRESHOLD 25
#endif

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

enum  BEH {//Behaviour
    LEFT, STAY, RIGHT, UNDEFINED, SLOW, FAST
  };

void next_points(int n, double mycar_s, int mycar_l, std::vector<double> &xpts, std::vector<double> &ypts, std::vector<double> &map_waypoints_s, std::vector<double> &map_waypoints_x, std::vector<double> &map_waypoints_y){
  for(int i=1;i<=n;i++){
    vector<double> next_pt = getXY(mycar_s+25*i, 4*(mycar_l)+2, map_waypoints_s, map_waypoints_x, map_waypoints_y);
    xpts.push_back(next_pt[0]);
    ypts.push_back(next_pt[1]);
  }
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  bool just_staring = true;
  double mycar_v_c = 0.0;
  int mycar_l = 1;

  h.onMessage([&mycar_l, &mycar_v_c, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double mycar_x = j[1]["x"];
          	double mycar_y = j[1]["y"];
          	double mycar_s = j[1]["s"];
          	double mycar_d = j[1]["d"];
          	double mycar_yaw = j[1]["yaw"];
          	double mycar_speed = j[1]["speed"];

            if(mycar_d > 0 && mycar_d < 4) mycar_l = 0;
            else if(mycar_d > 4 && mycar_d < 8) mycar_l = 1;
            else if(mycar_d > 8 && mycar_d < 12) mycar_l = 2;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];

          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

            int n = previous_path_x.size();
            if (n > 0) mycar_s = end_path_s;

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

            int car_l;
            bool car_lf = false, car_cf = false, car_rf = false, car_lb = false, car_cb = false, car_rb = false;
            for(int i=0; i < sensor_fusion.size(); i++){
              car_l = -1;
              double car_d = sensor_fusion[i][6];
              if(car_d > 0 && car_d < 4) car_l = 0;
              else if(car_d > 4 && car_d < 8) car_l = 1;
              else if(car_d > 8 && car_d < 12) car_l = 2;
              else{
                continue;
                std::cout << "Check for error" << std::endl;
              }
              double car_vx = sensor_fusion[i][3];
              double car_vy = sensor_fusion[i][4];
              double car_s = sensor_fusion[i][5];
              double car_speed = (car_vx*car_vx + car_vy*car_vy);
              car_speed += (double)n*TIME_INTERVAL*car_s;

              if (fabs(car_s - mycar_s) < DISTANCE_THRESHOLD){
                if (car_s >= mycar_s){
                  if (car_l-mycar_l==-1) car_lf=true;
                  else if (car_l-mycar_l==0) car_cf=true;
                  else if (car_l-mycar_l==1) car_rf=true;
                } //forward or back
                else{
                  if (car_l-mycar_l==-1) car_lb=true;
                  else if (car_l-mycar_l==0) car_cb=true;
                  else if (car_l-mycar_l==1) car_rb=true;
                }
              }
            }
            double delta_speed = 0;
            // Under the assumption that the lateral speed << forward speed
            if(car_cf){
              if(mycar_l > 0 && !car_lf && !car_lb) mycar_l-=1;
              else if(mycar_l < 2 && !car_rf && !car_rb) mycar_l+=1;
              else delta_speed -= MAX_ACC;
            }else{

#ifdef CONSERVATIVE
              // Prioritise center lane
              if((mycar_l != 1) && ((mycar_l == 0 && !car_rf && !car_rb) || (mycar_l == 2 && !car_lf && !car_lb))) mycar_l = 1; // Back to center.
#else
              // Prioritise inner lane
              if((mycar_l != 0) && (!car_lf && !car_lb)) mycar_l -=1; // Back to center.
#endif
              if (mycar_v_c < MAX_SPEED) delta_speed += MAX_ACC;
            }

            vector<double> xpts;
            vector<double> ypts;

            double start_x = mycar_x;
            double start_y = mycar_y;
            double start_yaw = deg2rad(mycar_yaw);
            double start_x_prev;
            double start_y_prev;
            if (n < 2){
                start_x_prev = start_x - cos(mycar_yaw);
                start_y_prev = start_y - sin(mycar_yaw);
            }else{ // Use last two points
                start_x = previous_path_x[n - 1];
                start_y = previous_path_y[n - 1];
                start_x_prev = previous_path_x[n - 2];
                start_y_prev = previous_path_y[n - 2];
                start_yaw = atan2(start_y-start_y_prev, start_x-start_x_prev);
            }
            xpts.push_back(start_x_prev);
            ypts.push_back(start_y_prev);
            xpts.push_back(start_x);
            ypts.push_back(start_y);

            next_points(N_POINTS_FOR_SPLINE, mycar_s, mycar_l, xpts, ypts, map_waypoints_s, map_waypoints_x, map_waypoints_y);

            for (int i = 0; i < xpts.size(); i++){
              double dx = xpts[i] - start_x;
              double dy = ypts[i] - start_y;
              xpts[i] = dx * cos(0 - start_yaw) - dy * sin(0 - start_yaw);
              ypts[i] = dx * sin(0 - start_yaw) + dy * cos(0 - start_yaw);
            }

            tk::spline sp;
            sp.set_points(xpts, ypts);

            for (int i = 0; i < n; i++){
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
            }

            double end_x = 25.0;
            double end_y = sp(end_x); // Get corresponding y point
            double end_dist = sqrt(end_x*end_x + end_y*end_y);

            double x_prev_main = 0;

            for(int i = 1; i < N_POINTS - n; i++){
              mycar_v_c += delta_speed; // Increment speed
              if (mycar_v_c > MAX_SPEED) mycar_v_c = MAX_SPEED;
              else if (mycar_v_c < MAX_ACC) mycar_v_c = MAX_ACC;

              double N = end_dist/(TIME_INTERVAL*mycar_v_c/CONVERT_M_2_MILES);
              double xpt = x_prev_main + end_x/N;
              double ypt = sp(xpt); // Get corresponding y point

              x_prev_main = xpt;

              double x_prev = xpt;
              double y_prev = ypt;

              xpt = start_x + x_prev * cos(start_yaw) - y_prev * sin(start_yaw);
              ypt = start_y + x_prev * sin(start_yaw) + y_prev * cos(start_yaw);
              next_x_vals.push_back(xpt);
              next_y_vals.push_back(ypt);
            }

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
