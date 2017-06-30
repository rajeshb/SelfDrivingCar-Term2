/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    // Number of particles
    num_particles = 50;
    
    default_random_engine gen;
    
    // Create normal distributions
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i = 0; i < num_particles; ++i) {
        Particle particle;
        particle.id = i;
        particle.x = dist_x(gen);
        particle.y = dist_y(gen);
        particle.theta = dist_theta(gen);
        particle.weight = 1.0;
        
        // Add to list of particles
        particles.push_back(particle);
        
        // Add to list of weights
        weights.push_back(particle.weight);
    }
    
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    
    default_random_engine gen;
    
    // Create normal distributions for noise
    normal_distribution<double> noise_x(0.0, std_pos[0]);
    normal_distribution<double> noise_y(0.0, std_pos[1]);
    normal_distribution<double> noise_theta(0.0, std_pos[2]);
    
    // Add measurement and random Gaussian noise to each of the particles
    for (auto && particle : particles) {
        // Avoid divide by zero
        if (yaw_rate == 0.0) {
            particle.x += velocity * delta_t * cos(particle.theta);
            particle.y += velocity * delta_t * sin(particle.theta);
            particle.theta = particle.theta;
        }
        else {
            particle.x += (velocity/yaw_rate) * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta));
            particle.y += (velocity/yaw_rate) * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t));
            particle.theta += yaw_rate * delta_t;
        }
        // Noise
        particle.x += noise_x(gen);
        particle.y += noise_y(gen);
        particle.theta += noise_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html
    
    // Prepare following calculations once, keep it ready
    double x2 = pow(std_landmark[0], 2);
    double y2 = pow(std_landmark[1], 2);
    double xy = std_landmark[0] * std_landmark[1];
    
    for (int i = 0; i < num_particles; ++i) {

        Particle& particle = particles[i];
        
        // weight is a product so init to 1.0
        double weight = 1.0;
        
        for (int j=0; j < observations.size(); j++) {
            
            // transform vehicle's observation to global coordinates
            LandmarkObs obs = observations[j];
            
            // predict landmark x, y. Equations from trigonometry.
            double predicted_x = obs.x * cos(particle.theta) - obs.y * sin(particle.theta) + particle.x;
            double predicted_y = obs.x * sin(particle.theta) + obs.y * cos(particle.theta) + particle.y;
            
            // initialise terms
            Map::single_landmark_s nearest_landmark;
            double min_distance = sensor_range;
            double distance = 0.0;
            
            // associate sensor measurements to map landmarks
            for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
                
                Map::single_landmark_s landmark = map_landmarks.landmark_list[k];
                
                // calculate distance between landmark and transformed observations
                distance = fabs(predicted_x - landmark.x_f) + fabs(predicted_y - landmark.y_f);
                
                // update nearest landmark to obs
                if (distance < min_distance) {
                    min_distance = distance;
                    nearest_landmark = landmark;
                }
                
            } // end associate nearest landmark
            
            // formula
            double x_diff = predicted_x - nearest_landmark.x_f;
            double y_diff = predicted_y - nearest_landmark.y_f;
            double numerator = exp( -0.5 * ( pow(x_diff, 2) / x2 + pow(y_diff, 2) / y2) );
            double denominator = 2.0 * M_PI * xy;
            
            // multiply particle weight by this obs-weight pair stat
            weight *= numerator / denominator;
            
        } // end observations loop
        
        // update particle weight
        particle.weight = weight;
        
        // update weight in PF array
        weights[i] = weight;
    }
}

void ParticleFilter::resample() {
    // Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    default_random_engine gen;
    
    // Take a discrete distribution equal to weights
    discrete_distribution<> disc_weights(weights.begin(), weights.end());
    
    // Resampled particles
    vector<Particle> resampled_particles;
    
    // Generate resampled particles
    for (int i = 0; i < num_particles; ++i)
        resampled_particles.push_back(particles[disc_weights(gen)]);
    
    // overwrite original with resampled particles
    particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    
    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
    
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
