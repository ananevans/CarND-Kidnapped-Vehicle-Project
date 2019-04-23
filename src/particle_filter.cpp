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


inline double multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
		double mu_x, double mu_y) {
	// calculate normalization term
	double gauss_norm;
	gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

	// calculate exponent
	double exponent;
	exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
            		   + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));

	// calculate weight using normalization terms and exponent
	double weight;
	weight = gauss_norm * exp(-exponent);

	return weight;
}


/**
 * Sets the number of particles.
 * Initializes all particles to first position (based on estimates of x, y, theta
 * and their uncertainties from GPS) and all weights to 1.
 * Adds random Gaussian noise to each particle.
 * Sets is_initialized to true.
 */
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	num_particles = 1000;  // TODO: Set the number of particles
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


/**
 * Adds measurements to each particle and add random Gaussian noise.
 */
void ParticleFilter::prediction(double delta_t, double std_pos[], 
		double velocity, double yaw_rate) {
	std::default_random_engine gen;
	for (int i = 0; i < num_particles; i++) {
		double x, y, theta;
		if ( yaw_rate == 0 ) {
			x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			theta = particles[i].theta;
		} else {
			x = particles[i].x + (velocity/yaw_rate) *
					(sin(particles[i].theta + delta_t*yaw_rate) - sin(particles[i].theta));
			y = particles[i].y + (velocity/yaw_rate) *
					(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			// error if I make theta between -pi, pi
			theta = particles[i].theta + yaw_rate * delta_t;
		}

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
	vector<double> new_weights;
	// update weights for each particle
	for ( int i = 0; i < num_particles; i++) { // for each particle
		//transform the observations from the vehicle to the world coordinate system wrt particle
		// transform to map x coordinate
		vector<LandmarkObs> transformed_obs;
		double x_part = particles[i].x;
		double y_part = particles[i].y;
		double theta = particles[i].theta;
		// create predicted vector of landmarks in sensor range
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
		if ( !observations.empty() ) {
			double w = 1.0; // new particle weight
			vector<int> associations;
			vector<double> sense_x;
			vector<double> sense_y;
			for ( auto obs = observations.begin(); obs != observations.end(); obs++ ) {
				// transform the observation coordinates into the particle's frame
				double x_obs = obs->x;
				double y_obs = obs->y;
				LandmarkObs tobs;
				tobs.x = x_part + (cos(theta) * x_obs) - (sin(theta) * y_obs);
				tobs.y = y_part + (sin(theta) * x_obs) + (cos(theta) * y_obs);
				// find the associated landmark to this observation
				LandmarkObs* closest_landmark = NULL;
				double min_dist = sensor_range;
				for ( auto pred = predicted.begin(); pred != predicted.end(); pred++) {
					double d = dist(tobs.x, tobs.y, pred->x, pred->y);
					if (min_dist > d) {
						closest_landmark = &(*pred);
						min_dist = d;
					}
				}
				if ( closest_landmark != NULL ) {
					// update weight
					double w1 = multiv_prob(std_landmark[0], std_landmark[1],
							tobs.x, tobs.y,
							closest_landmark->x, closest_landmark->y);
					if ( w1 != 0 ) {
						w *= w1;
					}
					// update associations for this particle
					associations.push_back( closest_landmark->id );
				} else {
					std::cout << "Can not find landmark to associate to observation ("
							<< obs->x << ","
							<< obs->y << ")";
					associations.push_back(1);
				}
				// update the coordinates for this association
				sense_x.push_back(tobs.x);
				sense_y.push_back(tobs.y);
			}
			particles[i].weight = w;
			new_weights.push_back(w);
			SetAssociations(particles[i], associations, sense_x, sense_y);
		} else {
			std::cout << "NO OBSERVATIONS!" << std::endl;
			new_weights.push_back(particles[i].weight);
		}

	}
	weights = new_weights;
}

/**
 * Resample particles with replacement with probability proportional
 *   to their weight.
 */
void ParticleFilter::resample() {
	std::vector<Particle> new_particles;
	std::default_random_engine gen;
	std::discrete_distribution<int> d(weights.begin(), weights.end());
	for ( int i = 0; i < num_particles; i++) {
		new_particles.push_back(particles[d(gen)]);
	}
	particles = new_particles;
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
