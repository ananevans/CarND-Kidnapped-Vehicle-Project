/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 5000;  // TODO: Set the number of particles
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  for (int i = 0; i < num_particles; i++) {
	  Particle p;
	  p.id = i;
	  p.x = dist_x(gen);
	  p.y = dist_y(gen);
	  p.theta = dist_theta(gen);
	  p.weight = 1.0;
	  particles.push_back(p);
	  weights.push_back(p.weight);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
	std::default_random_engine gen;
	for (int i = 0; i < num_particles; i++) {
		double x = particles[i].x + (velocity/yaw_rate) *
				(sin(particles[i].theta + delta_t*yaw_rate) - sin(particles[i].theta));
		double y = particles[i].y + (velocity/yaw_rate) *
				(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
		// error if I make theta between -pi, pi
		double theta = particles[i].theta + yaw_rate * delta_t;
		normal_distribution<double> dist_x(x, std_pos[0]);
		normal_distribution<double> dist_y(y, std_pos[1]);
		normal_distribution<double> dist_theta(theta, std_pos[2]);

		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
	// assumption: observations and predicted are in the same coordinate system
	for ( auto obs = observations.begin(); obs != observations.end(); obs++ ) {
		// Find closest prediction
		int id = -1;
		double min_dist = std::numeric_limits<double>::max();
		for ( auto pred = predicted.begin(); pred != predicted.end(); pred++) {
			double d = dist(obs->x, obs->y, pred->x, pred->y);
			if (min_dist > d) {
				id = pred->id;
				min_dist = d;
			}
		}
		if ( id >= 0 ) {
			obs->id = id;
		} else {
			obs->id = -1;
			std::cout << "Prediction not found for observation "
					<< obs->x << " , " << obs->y
					<<"!" << std::endl;
		}
	}
}

/**
 * updateWeights Updates the weights for each particle based on the likelihood
 *   of the observed measurements.
 * @param sensor_range Range [m] of sensor
 * @param std_landmark[] Array of dimension 2
 *   [Landmark measurement uncertainty [x [m], y [m]]]
 * @param observations Vector of landmark observations
 * @param map Map class containing map landmarks
 */
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

	vector<double> new_weights;
	// update weights for each particle
	for ( int i = 0; i < num_particles; i++) { // for each particle
		//transform the observations from the vehicle to the world coordinate system wrt particle
		// transform to map x coordinate
		vector<LandmarkObs> transformed_obs;
		double x_part = particles[i].x;
		double y_part = particles[i].y;
		double theta = particles[i].theta;
		for ( auto obs = observations.begin(); obs != observations.end(); obs++ ) {
			double x_obs = obs->x;
			double y_obs = obs->y;
			LandmarkObs tobs;
			tobs.x = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
			tobs.y = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
			tobs.id = obs->id;
			transformed_obs.push_back(tobs);
		}
		// create predicted vector
		vector<LandmarkObs> predicted;
		for ( auto ml = map_landmarks.landmark_list.begin(); ml != map_landmarks.landmark_list.end(); ml++ ) {
			if ( dist(ml->x_f, ml->y_f, particles[i].x, particles[i].y) < sensor_range ) {
				LandmarkObs p;
				p.x = ml->x_f;
				p.y = ml->y_f;
				p.id = ml->id_i;
				predicted.push_back(p);
			}
		}
		// associate measurements with landmarks
		dataAssociation( predicted, transformed_obs);
		//determine measurement probability
		if ( !transformed_obs.empty() ) {
			double w = 1.0;
			for ( auto tobs = transformed_obs.begin(); tobs != transformed_obs.end(); tobs++ ) {
				// find the associated landmark
				for ( auto pred = predicted.begin(); pred != predicted.end(); pred++ ) {
					if (tobs->id == pred->id) {
						w *= multiv_prob(std_landmark[0], std_landmark[1],
												tobs->x, tobs->y,
												pred->x, pred->y);
						break;
					}
				}
			}
			particles[i].weight = w;
			new_weights.push_back(w);
		} else {
			std::cout << "NO OBSERVATIONS!" << std::endl;
			new_weights.push_back(particles[i].weight);
		}
	}
	weights = new_weights;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
	std::vector<Particle> p1;
	std::vector<double> w1;
	std::default_random_engine gen;
	std::discrete_distribution<int> d(0,num_particles);

	int index = d(gen);
	double w = *std::max_element(weights.begin(), weights.end());

	std::uniform_real_distribution<double> u(0, 2*w);

	double beta = 0.0;
	for ( int i = 0; i < num_particles; i++) {
		beta = beta + u(gen);
		while ( weights[index] < beta ) {
			beta = beta - weights[index];
			index = (index+1) % num_particles;
		}
		p1.push_back(particles[index]);
		w1.push_back(particles[index].weight);
	}

	particles = p1;
	weights = w1;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
