#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <list>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  cout << "map_waypoints_x SIZE: " << map_waypoints_x.size() << endl;
  
  
  /*
  vector<double> mw_x;
  vector<double> mw_y;
  vector<double> mw_s;

  for (int i = 1; i < (map_waypoints_x.size() - 5); i++) {

	 // DETERMINE POINTS OF THE SPLINE 

	  vector<double> mw_x_spline;
	  vector<double> mw_y_spline;
	  double ref_x_wp;
	  double ref_y_wp;
	  double ref_yaw_wp;
	  vector<double> ptsx_wp;
	  vector<double> ptsy_wp;
	  double ptsx_wp_future;
	  double ptsy_wp_future;
	  
	  ref_x_wp = map_waypoints_x[i];
	  ref_y_wp = map_waypoints_y[i];
	  ref_yaw_wp = atan2(ref_y_wp - map_waypoints_y[i-1], ref_x_wp - map_waypoints_x[i-1]);

	  for (int j = 0; j < 3; j++) {
		  // Shift car reference angle to 0 degrees
		  double shift_x_wp = map_waypoints_x[i+j+1] - ref_x_wp;
		  double shift_y_wp = map_waypoints_y[i+j+1] - ref_y_wp;
		  
		  ptsx_wp_future = (shift_x_wp * cos(0 - ref_yaw_wp) - shift_y_wp * sin(0 - ref_yaw_wp));
		  ptsy_wp_future = (shift_x_wp * sin(0 - ref_yaw_wp) + shift_y_wp * cos(0 - ref_yaw_wp));

		  ptsx_wp.push_back(ptsx_wp_future);
		  ptsy_wp.push_back(ptsy_wp_future);
		  		 
		  mw_x_spline.push_back(ptsx_wp_future);
		  mw_y_spline.push_back(ptsy_wp_future);
	  }
	 
	  // Create a spline called map_wp
	  tk::spline map_wp;

	  // Set (x,y) points to the spline
	  map_wp.set_points(mw_x_spline, mw_y_spline);

	  // DETERMINE NEW WAYPOINTS FROM SPLINE

	  // Compute how to break up spline points so we travel at our desired reference velocity
	  double target_x = ptsx_wp[1];
	  double target_y = ptsy_wp[1];
	  double target_dist = sqrt((target_x - ptsx_wp[0]) * (target_x - ptsx_wp[0]) + (target_y - ptsy_wp[0]) * (target_y - ptsy_wp[0])); 
	  double x_add_on = 0;
	  double target_s = map_waypoints_s[i + 1];
	  double target_dist_s = map_waypoints_s[i + 1] - map_waypoints_s[i];
	  double s_add_on = 0;

	  // Fill up the rest of the path planner to always output 30 points
	  for (int k = 0; k < 30; k++) {
		  double N = 30;
		  double x_point_car = ptsx_wp[0] + x_add_on + (target_x - ptsx_wp[0]) / N;
		  double y_point_car = map_wp(x_point_car);
		  //cout << "ptsx_wp[0]: " << ptsx_wp[0] << endl;
		  x_add_on = x_point_car - ptsx_wp[0];

		  // Come back to global coordinates
		  double x_ref_wp = x_point_car;
		  double y_ref_wp = y_point_car;

		  // Rotate back to normal after rotating it earlier
		  double x_point = (x_ref_wp * cos(ref_yaw_wp) - y_ref_wp * sin(ref_yaw_wp));
		  double y_point = (x_ref_wp * sin(ref_yaw_wp) + y_ref_wp * cos(ref_yaw_wp));

		  x_point += ref_x_wp;
		  y_point += ref_y_wp;

		  mw_x.push_back(x_point);
		  mw_y.push_back(y_point);

		  double s_point = map_waypoints_s[i] + s_add_on + (map_waypoints_s[i + 1] - map_waypoints_s[i]) / N;
		  s_add_on = s_point - map_waypoints_s[i];

		  mw_s.push_back(s_point);
		  
	  }
  }
  */

  // start in lane 1
  int lane = 1;

  int previous_lane = 1;

  //target velocity
  double ref_vel = 0.0;

  bool lane_changing = false;

  bool move_two_lanes = false;

  int time_change = 0;

  int two_lane_change_starting_lane = 0;
    

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy, &lane, &previous_lane, &ref_vel, &lane_changing, &move_two_lanes, &time_change, &two_lane_change_starting_lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
			vector<double> test;


          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

			int prev_size = previous_path_x.size();

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

			// CHECKING OTHER CARS

			if (prev_size >0 ) {
				car_s = end_path_s;
			}

			bool too_close = false;
			vector<bool> lane_free(3, true);
			vector<bool> lane_available(3, true);
			
			vector<float> mylist;

			//cout << "D position of other cars: " << endl;
			// find ref_v to use
			for (int i=0; i < sensor_fusion.size(); i++) {
				// car is in my lane
				double d = sensor_fusion[i][6];
				double check_car_s = sensor_fusion[i][5];
				
				//vector<double> pos_car(2);
				//pos_car.push_back(check_car_s);
				//mylist.push_back(d);
				//vector<double> cars;
				//cars.push_back(pos_car);
						
				double vx = sensor_fusion[i][3];
				double vy = sensor_fusion[i][4];
				double check_speed = sqrt(vx*vx + vy * vy);


				check_car_s += ((double)prev_size * .02 * check_speed); //  if using previous points can project a value outward
																		
				// consider Average Dry Stopping Distance (Assuming two-thirds of a second reaction time)
				// https://en.wikibooks.org/wiki/Driving/Safety/Wet_weather_driving

				double safety_dist = (ref_vel + ((ref_vel) * (ref_vel) / 20)) * 0.3048;

				// CHECK COMPLETELY FREE LINES
				if ((check_car_s > car_s - 25) && (check_car_s - car_s < 100)) {
					if ((d < 4) && (d > 0)) {
						lane_free[0] = false;
					}
					else if ((d >= 4) && (d < 8)) {
						lane_free[1] = false;
					}
					else if ((d >= 8) && (d < 12)) {
						lane_free[2] = false;
					}
				}

				// CHECK AVAILABLE LINES
				if ((check_car_s > (car_s - 10)) && (check_car_s - car_s < safety_dist)) {

					
					if ((d < 4) && (d > 0)) {
						lane_available[0] = false;
					}
					else if ((d >= 4) && (d < 8)) {
						lane_available[1] = false;
					}
					else if ((d >= 8) && (d < 12)) {
						lane_available[2] = false;
					}
				}

				// check s values greater than mine and s gap
				if (d < (2 + 4 * lane + 3) && d >(2 + 4 * lane - 3)) {
					if ((check_car_s > car_s) && (check_car_s - car_s < safety_dist)) {

						// DO SOMETHING TO AVOID COLLISION
						too_close = true;
					}
				}

			}
			
			/*
			// FREE LINES

			cout << "FREE LINES: ";

			if (lane_free[0] == true) {
				cout << "LANE LEFT, ";
			}
			if (lane_free[1] == true) {
				cout << "LANE CENTER, ";
			}
			if (lane_free[2] == true) {
				cout << "LANE RIGHT, ";
			}
			cout << endl;
			*/

			// AVAILABLE LINES

			/*
			cout << "AVAILABLE LANES: ";

			if (lane_available[0] == true) {
				cout << "LANE LEFT, ";
			}
			if (lane_available[1] == true) {
				cout << "LANE CENTER, ";
			}
			if (lane_available[2] == true) {
				cout << "LANE RIGHT, ";
			}
			cout << endl;
			*/

			// DECREASE SPEED

			if (too_close) {
				ref_vel -= .25;
			} else if ((move_two_lanes) && ref_vel > 30) {
				ref_vel -= .04;
			}
			else if (ref_vel < 49.75) {
				ref_vel += .25;
			}
			
			// CONSIDER TIME FOR LANE CHANGING
			if (lane_changing) {
				time_change = time_change + 1;
				move_two_lanes = false;
				cout << "time_change: " << time_change << endl;
			}

			if (time_change > 150) {
				lane_changing = false;
				time_change = 0;
			}

			// CONSIDER TIME FOR TWO LANE CHANGING SITUATION
			if (move_two_lanes) {
				cout << "Trying to move to the opposite lane" << endl;
			}

			previous_lane = lane;

			// LANE CHANGING POSSIBILE SITUATIONS
			// If the center lane is considered free, move to the center lane
			if (lane_free[1]) {
				lane = 1;
				if (previous_lane != lane) {
					lane_changing = true;
				}			
			} else if (((too_close) && (lane_changing == false)) || (move_two_lanes)){
				lane_changing = true;
				if (lane == 0 ) {
					if (lane_available[1] == true) {
						lane++ ;
					} else  {
						if (lane_available[2] == true) {
							// Should I move to the right lane without touching any other car
							move_two_lanes = true;
							lane_changing = false;
						} else {
							cout << "Wait..." << endl;
							lane_changing = false;
							move_two_lanes = false;
						}
					}
				} else if (lane == 1) {
					if ((lane_free[2]) && (lane_free[0] = false)) {
						lane++ ;
					} else if (lane_available[0] == true) {
						lane-- ;
					} else if (lane_available[2] == true) {
						lane++ ;
					} else {
						cout << "Wait..." << endl;
						lane_changing = false;
					}

				} else if (lane == 2) {
					if (lane_available[1] == true) {
						lane-- ;
					} else if (lane_available[1] == false) {
						if (lane_available[0] == true) {
							// Should I move to the left lane without touching any other car
							move_two_lanes = true;
							lane_changing = false;
						} else {
							cout << "Wait..." << endl;
							lane_changing = false;
							move_two_lanes = false;
						}
					}
				} else {
					cout << "ERROR" << endl;
				}
			}

			// Reset to initiale state move_two_lanes variable
			if (previous_lane != lane) {
			move_two_lanes = false;
			}

			if (lane_free[lane]) {
				move_two_lanes = false;
			}

			//END

          	json msgJson;

			// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
			double next_d_vals;
			double next_s_vals;

			vector<double> ptsx;
			vector<double> ptsy;

			double ref_x = car_x;
			double ref_y = car_y;
			double ref_yaw = deg2rad(car_yaw);

			// if previous size is almost empty, use the car as a starting reference
			if (prev_size < 2) {
				// use two points that make the path tangent to the car
				double prev_car_x = car_x - cos(car_yaw);
				double prev_car_y = car_y - sin(car_yaw);

				ptsx.push_back(prev_car_x);
				ptsx.push_back(car_x);
				
				ptsy.push_back(prev_car_y);
				ptsy.push_back(car_y);
			}

			// use the prevoius path's end point as starting reference
			else {
				// redefine reference state as previous path and point
				ref_x = previous_path_x[prev_size - 1];
				ref_y = previous_path_y[prev_size - 1];

				double ref_x_prev = previous_path_x[prev_size - 2];
				double ref_y_prev = previous_path_y[prev_size - 2];
				ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);
				
				//use two points that make the path tangent to the previous path's end point
				ptsx.push_back(ref_x_prev);
				ptsx.push_back(ref_x);

				ptsy.push_back(ref_y_prev);
				ptsy.push_back(ref_y);
			}

			// Using Frenet, add 30 m evenly spaced points ahead of the starting reference
			vector<double> next_wp0 = getXY(car_s + 30, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp1 = getXY(car_s + 60, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			vector<double> next_wp2 = getXY(car_s + 90, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

			ptsx.push_back(next_wp0[0]);
			ptsx.push_back(next_wp1[0]);
			ptsx.push_back(next_wp2[0]);

			ptsy.push_back(next_wp0[1]);
			ptsy.push_back(next_wp1[1]);
			ptsy.push_back(next_wp2[1]);

			for (int i = 0; i < ptsx.size(); i++) {
				// Shift car reference angle to 0 degrees
				double shift_x = ptsx[i] - ref_x;
				double shift_y = ptsy[i] - ref_y;

				ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
				ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
			}


			// Create a spline called s
			tk::spline s;

			// Set (x,y) points to the spline
			s.set_points(ptsx, ptsy);

			// Start with all the previous path points from last time
			for (int i = 0; i < previous_path_x.size(); i++) {
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}

			// Compute how to break up spline points so we travel at our desired reference velocity
			double target_x = 30.0;
			double target_y = s(target_x);
			double target_dist = sqrt((target_x) * (target_x)+(target_y) * (target_y));
			double x_add_on = 0;

			// Fill up the rest of the path planner to always output 50 points
			for (int i = 1; i <= 50 - previous_path_x.size(); i++) {
				double N = (target_dist / (.02*ref_vel / 2.24));
				double x_point = x_add_on + (target_x) / N;
				double y_point = s(x_point);

				x_add_on = x_point;

				double x_ref = x_point;
				double y_ref = y_point;

				// Rotate back to normal after rotating it earlier
				x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
				y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

				x_point += ref_x;
				y_point += ref_y;

				next_x_vals.push_back(x_point);
				next_y_vals.push_back(y_point);
			}
			
			// END

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
