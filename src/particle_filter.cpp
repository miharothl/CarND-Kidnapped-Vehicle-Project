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
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  normal_distribution<double> N_x(x, std[0]);
  normal_distribution<double> N_y(y, std[1]);
  normal_distribution<double> N_theta(theta, std[2]);

    default_random_engine gen;

  this->num_particles = 30;

  for (auto i=0; i<num_particles; i++)
  {
    auto particle = Particle();
    particle.id = i;
    particle.x = N_x(gen);
    particle.y = N_y(gen);
    particle.theta = N_theta(gen);
    particle.weight = 1;

    this->particles.push_back(particle);
  }

  this->is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/


  normal_distribution<double> N_x(0, std_pos[0]);
  normal_distribution<double> N_y(0, std_pos[1]);
  normal_distribution<double> N_theta(0, std_pos[2]);

  if (fabs(yaw_rate) > 0.00001)
  {
    for (auto &particle : particles ) {
      particle.x = particle.x + (velocity / yaw_rate) * (sin(particle.theta + yaw_rate * delta_t) - sin(particle.theta));
      particle.y = particle.y + (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + yaw_rate * delta_t));
      particle.theta = particle.theta + (yaw_rate*delta_t);

      particle.x += N_x(gen);
      particle.y += N_y(gen);
      particle.theta += N_theta(gen);
    }
  }
  else
  {
    for (auto &particle : particles ) {
      particle.x = particle.x + (velocity * delta_t * cos(particle.theta));
      particle.y = particle.y + (velocity * delta_t * sin(particle.theta));

      particle.x += N_x(gen);
      particle.y += N_y(gen);
      particle.theta += N_theta(gen);
    }
  }
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.

  for (LandmarkObs& p : observations) {


    float min_d = 10000;

    for (LandmarkObs& r : predicted) {
      float d = dist(p.x, p.y, r.x, r.y);

      if (d < min_d) {
        min_d = d;
        p.id = r.id;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
    std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation
  //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account
  //   for the fact that the map's y-axis actually points downwards.)
  //   http://planning.cs.uiuc.edu/node99.html

  vector<LandmarkObs> landmarks;

  for (Map::single_landmark_s sl : map_landmarks.landmark_list) {

    LandmarkObs l;
    l.id = sl.id_i;
    l.x = sl.x_f;
    l.y	 = sl.y_f;

    landmarks.push_back(l);
  }



  this->weights.clear();

  for (Particle& particle : this->particles) {

    vector<LandmarkObs> transformedObservations;

    for (LandmarkObs observation : observations) {

//			float tranObsX = particle.x + cos(particle.theta) * observation.x + sin(particle.theta) * observation.y;
//			float tranObsY = particle.y - sin(particle.theta) * observation.x + cos(particle.theta) * observation.y;


      float tranObsX = particle.x + cos(particle.theta) * observation.x - sin(particle.theta) * observation.y;
      float tranObsY = particle.y + sin(particle.theta) * observation.x + cos(particle.theta) * observation.y;

          LandmarkObs tranObs;
      tranObs.id = observation.id;
          tranObs.x = tranObsX;
      tranObs.y = tranObsY;

      transformedObservations.push_back(tranObs);
    }

    this->dataAssociation(landmarks, transformedObservations);



    vector<int> assoc;
    vector<double> sense_x;
    vector<double> sense_y;


    vector<double> p_xy;
    for (LandmarkObs& to : transformedObservations) {
      to.id;

      LandmarkObs found;
      for (LandmarkObs l : landmarks) {

        if (l.id == to.id) {
          found = l;
          assoc.push_back(found.id);
          sense_x.push_back(found.x);
          sense_y.push_back(found.y);

        }
      }

      double a= pow(to.x - found.x, 2) / (2 * pow(std_landmark[0], 2));
      double b= pow(to.y - found.y, 2) / (2 * pow(std_landmark[1], 2));

      double ex = exp(-( a+b ));

      double p = 1 / (2*3.14*std_landmark[0]*std_landmark[1]) * ex;

      p_xy.push_back(p);


    }

//		this->SetAssociations(particle, assoc, sense_x, sense_y);

    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations= assoc;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    particle.weight = 1;

    for (double p : p_xy) {
      particle.weight = particle.weight * p;
    }

    this->weights.push_back(particle.weight);

    if (particle.weight != 0) {
      int break_me = 0;

      std::cout << "weight:" << particle.id << ": = " << particle.weight << endl;
    }
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight.
  // NOTE: You may find std::discrete_distribution helpful here.
  // http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::discrete_distribution<int> D_weights(this->weights.begin(), this->weights.end());

  std::vector<Particle> resampled_particles;

  for (auto particle : this->particles) {
    int id = D_weights(gen);
    resampled_particles.push_back(this->particles[id]);
  }

  this->particles = resampled_particles;
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
