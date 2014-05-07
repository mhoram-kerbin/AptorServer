/**
 * Launch.cpp - A Launch Trajectory Optimizer for KSP based on PSOPT
 *
 * This file is part of AptorServer.
 *
 * Copyright 2014 Mhoram Kerbin
 *
 * AptorServer is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License,
 * or (at your option) any later version.
 *
 * AptorServer is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AptorServer. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "Launch.h"

#include "Exception/ServerException.h"

#include <psopt.h>
/**
 * Assumptions:
 *
 * Time always starts at 0.
 *
 */

namespace Aptor {


using std::vector;

#define NR_OF_STATES 7
#define NR_OF_CONTROLS 3


Launch::Launch(): Optimizer()
{
	commandMapper["NAME"]                   = [this] (vector<string> c) -> string { setName(c);                   return "OK"; };
	commandMapper["PLANET_MASS"]            = [this] (vector<string> c) -> string { setPlanetMass(c);             return "OK"; };
	commandMapper["PLANET_RADIUS"]          = [this] (vector<string> c) -> string { setPlanetRadius(c);           return "OK"; };
	commandMapper["PLANET_SCALE_HEIGHT"]    = [this] (vector<string> c) -> string { setPlanetScaleHeight(c);      return "OK"; };
	commandMapper["PLANET_P0"]              = [this] (vector<string> c) -> string { setPlanetP0(c);               return "OK"; };
	commandMapper["PLANET_ROTATION_PERIOD"] = [this] (vector<string> c) -> string { setPlanetRotationPeriod(c);   return "OK"; };
	commandMapper["PLANET_SOI"]             = [this] (vector<string> c) -> string { setPlanetSOI(c);              return "OK"; };
	commandMapper["LAUNCH_LATITUDE"]        = [this] (vector<string> c) -> string { setLaunchLatitude(c);         return "OK"; };
	commandMapper["LAUNCH_LONGITUDE"]       = [this] (vector<string> c) -> string { setLaunchLongitude(c);        return "OK"; };
	commandMapper["LAUNCH_ALTITUDE"]        = [this] (vector<string> c) -> string { setLaunchAltitude(c);         return "OK"; };
	commandMapper["MAX_VELOCITY"]           = [this] (vector<string> c) -> string { setShipMaxV(c);               return "OK"; };
	commandMapper["ADD_STAGE"]              = [this] (vector<string> c) -> string { addStage(c);                  return "OK"; };
	commandMapper["ADD_ENGINE"]             = [this] (vector<string> c) -> string { addEngine(c);                 return "OK"; };
	commandMapper["TARGET_PERIAPSIS"]       = [this] (vector<string> c) -> string { setTargetPeriapsis(c);        return "OK"; };
	commandMapper["COMPUTE"]                = [this] (vector<string> c) -> string { compute(c);                   return "OK"; };
	commandMapper["POSTPROCESS"]            = [this] (vector<string> c) -> string { performPostprocessing(true);  return "OK"; };
	commandMapper["VERBOSE_POSTPROCESS"]    = [this] (vector<string> c) -> string { performPostprocessing(false); return "OK"; };
	commandMapper["TEST"]                   = [this] (vector<string> c) -> string { Test();                       return "OK"; };
	commandMapper["GET_PITCH_THRUST"] = [this] (vector<string> c) -> string { return getPitchThrust(c); };
	commandMapper["GET_FINAL_TIMES"]  = [this] (vector<string> c) -> string { return getFinalPhaseTimes(c); };
	setupDaeFunction();
}


void Launch::setName(vector<string> command)
{
	configuration.name = command[1];
}
void Launch::setPlanetMass(vector<string> command)
{
	configuration.planetMass = getArgAsDouble(1, command);
	configuration.planetMu = configuration.planetMass * constants.gravitationalConstant;
}

void Launch::setPlanetRadius(vector<string> command)
{
	configuration.planetRadius = getArgAsDouble(1, command);
}
void Launch::setPlanetScaleHeight(vector<string> command)
{
	configuration.planetScaleHeight = getArgAsDouble(1, command);
}
void Launch::setPlanetP0(vector<string> command)
{
	configuration.planetP0 = getArgAsDouble(1, command);
}
void Launch::setPlanetRotationPeriod(vector<string> command)
{
	configuration.planetRotationPeriod = getArgAsDouble(1, command);
}
void Launch::setPlanetSOI(vector<string> command)
{
	configuration.planetSOI = getArgAsDouble(1, command);
}
void Launch::setLaunchLatitude(vector<string> command)
{
	configuration.launchLatitude = getArgAsDouble(1, command);
}
void Launch::setLaunchLongitude(vector<string> command)
{
	configuration.launchLongitude = getArgAsDouble(1, command);
}
void Launch::setLaunchAltitude(vector<string> command)
{
	configuration.launchAltitude = getArgAsDouble(1, command);
}
void Launch::setShipMaxV(vector<string> command)
{
	configuration.shipMaxV = getArgAsDouble(1, command);
}
void Launch::addStage(vector<string> command)
{
	StageParameters para;
	para.mass    = getArgAsDouble(1, command) * 1000;
	para.fuel    = getArgAsDouble(2, command) * 1000;
	para.drag    = getArgAsDouble(3, command);
	configuration.stages.push_back(para);
}
void Launch::addEngine(vector<string> command)
{
	double thrust = getArgAsDouble(1, command);
	double isp0 = getArgAsDouble(2, command);
	double ispv = getArgAsDouble(3, command);
	int nr = 1;
	if (command.size() > 4) {
		nr = getArgAsInt(4, command);
	}
	if (thrust <= 0 || isp0 <= 0 || ispv <= 0 || nr <= 0) {
		throw ServerException("Input arguments for ADD_ENGINE must be positive");
	}
	thrust = thrust * nr;
	StageParameters sp = configuration.stages.back();
	configuration.stages.back().combinedThrust += thrust;
	int i = 0;
	while (i < (int)configuration.stages.back().engines.size()) {
		if (configuration.stages.back().engines[i].isp1Atm == isp0 &&
				configuration.stages.back().engines[i].ispVac == ispv) {
			configuration.stages.back().engines[i].thrust += thrust;
			break;
		}
		i++;
	}
	if (i == (int)configuration.stages.back().engines.size()) {
		Engine e(thrust, isp0, ispv);
		configuration.stages.back().engines.push_back(e);
	}
}
void Launch::setTargetPeriapsis(vector<string> command)
{
	configuration.targetPeriapsis = getArgAsDouble(1, command);
}

void Launch::compute(vector<string> command)
{
	setupInitialState();
	setupTime();
	setupDae();
	setupEvents();
	setupLinkages();
	setupCosts();
	Optimizer::compute(configuration.name, configuration.name, calcLinkagesNumber());
	initResults();
}

void Launch::setupInitialState()
{
	double dist = configuration.planetRadius + configuration.launchAltitude;
	double f = dist * 2 * pi * cos(configuration.launchLatitude) / configuration.planetRotationPeriod;
	initialState.posX = dist * cos(configuration.launchLatitude) * cos(configuration.launchLongitude);
	initialState.posY = dist * cos(configuration.launchLatitude) * sin(configuration.launchLongitude);
	initialState.posZ = dist * sin(configuration.launchLatitude);
	initialState.velX = -sin(configuration.launchLongitude) * f;
	initialState.velY =  cos(configuration.launchLongitude) * f;
	initialState.velZ = 0;
	initialState.mass = configuration.stages[0].mass;
}

void Launch::setupTime()
{
	int stages = configuration.stages.size();
	double pressure = calcPressure(configuration.launchAltitude).value();
	timLower.push_back(0);
	timUpper.push_back(0);
	for (int stage = 1; stage <= stages; stage++) {
		double thrust = configuration.stages[stage-1].combinedThrust;
//		double isp = calcIsp(pressure, configuration.stages[stage-1].isp1Atm, configuration.stages[stage-1].ispVac).value();
//		double isp = calcIsp(pressure, configuration.stages[stage-1].engines).value();
		double isp = calcIsp(pressure, stage-1).value();
//		cout << "fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\nstage " << stage << " isp " << isp << " thrust " << thrust << endl;
		double fuelConsumption = thrust * 1000 / (isp * constants.g0); // conversion kN -> N
		timLower.push_back(timLower[stage-1] + configuration.stages[stage-1].fuel / fuelConsumption);
		timUpper.push_back(timUpper[stage-1] + configuration.stages[stage-1].fuel / fuelConsumption * (1 + 0.3 * stage));
		if (stage < stages) {
			psoptInterface.SetTimeConstraints(stage,
					timLower[stage-1] * 0.99, timUpper[stage-1],
					timLower[stage] * 0.99, timUpper[stage]);
		} else {
			double upper = 2 * std::max(timUpper[stage], 2 * configuration.timeToOrbit);
			psoptInterface.SetTimeConstraints(stage,
					timLower[stage-1] * 0.99, timUpper[stage-1],
					timLower[stage] * 0.99, upper);
		}
	}
}

void Launch::setupDaeFunction()
{
	daeFunction = [this] (adouble* derivatives, adouble* path, adouble* states,
			adouble* controls, adouble* parameters, adouble& time,
			adouble* xad, int iphase) -> void
			{
		adouble pos[3], vel[3];
		pos[0] = states[stateIndex[iphase-1].posX];
		pos[1] = states[stateIndex[iphase-1].posY];
		pos[2] = states[stateIndex[iphase-1].posZ];
		vel[0] = states[stateIndex[iphase-1].velX];
		vel[1] = states[stateIndex[iphase-1].velY];
		vel[2] = states[stateIndex[iphase-1].velZ];
		adouble pos2 = dot(pos, pos, 3);
		adouble distance = sqrt(pos2);
		double thrustFactor = 1000;
		adouble pressure = calcPressure(distance - configuration.planetRadius);
		adouble gvv[3];
		calcGroundVelocityVector(pos, vel, gvv);
		adouble density = pressure * constants.conversionFactor;
		adouble dragFactor = - 0.5 * density * configuration.stages[iphase-1].drag * 0.008 * states[stateIndex[iphase-1].mass] * sqrt(dot (gvv,gvv,3));
		adouble gravityFactor = - states[stateIndex[iphase-1].mass] * configuration.planetMu / (pos2 * distance);
		adouble thrust[3] = {controls[daeControlIndex[iphase-1].thrustX], controls[daeControlIndex[iphase-1].thrustY], controls[daeControlIndex[iphase-1].thrustZ]};
		adouble thrust2 = dot(thrust, thrust, 3);

		daeAccelX = (thrust[0] * thrustFactor + gvv[0] * dragFactor + pos[0] * gravityFactor) / states[stateIndex[iphase-1].mass];
		daeAccelY = (thrust[1] * thrustFactor + gvv[1] * dragFactor + pos[1] * gravityFactor) / states[stateIndex[iphase-1].mass];
		daeAccelZ = (thrust[2] * thrustFactor + gvv[2] * dragFactor + pos[2] * gravityFactor) / states[stateIndex[iphase-1].mass];
//		adouble isp = calcIsp(pressure, configuration.stages[iphase-1].isp1Atm, configuration.stages[iphase-1].ispVac);
//		adouble isp = calcIsp(pressure, configuration.stages[iphase-1].engines);
		adouble isp = calcIsp(pressure, iphase-1);
		derivatives[stateIndex[iphase-1].mass] = -sqrt(thrust2) / (isp * constants.g0) * 1000;

		derivatives[stateIndex[iphase-1].posX] = vel[0];
		derivatives[stateIndex[iphase-1].posY] = vel[1];
		derivatives[stateIndex[iphase-1].posZ] = vel[2];
		derivatives[stateIndex[iphase-1].velX] = daeAccelX;
		derivatives[stateIndex[iphase-1].velY] = daeAccelY;
		derivatives[stateIndex[iphase-1].velZ] = daeAccelZ;

		path[daePathIndex.thrust2] = thrust2;
		path[daePathIndex.distance2] = pos2;

		adouble ev[3];
		calcEccentricityVector(pos, vel, ev);
		path[daePathIndex.eccentricity2] = dot(ev, ev, 3);
		//printf("(%f,%f,%f) (%f,%f,%f) %f (%f,%f,%f)", states[0].value(), states[1].value(), states[2].value(), states[3].value(), states[4].value(), states[5].value(), states[6].value(), controls[0].value(), controls[1].value(), controls[2].value());
		//printf("(%f,%f,%f) (%f,%f,%f) %f\n", derivatives[0].value(), derivatives[1].value(), derivatives[2].value(), derivatives[3].value(), derivatives[4].value(), derivatives[5].value(), derivatives[6].value());
/*
		cout << "Pos (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << endl;
		cout << "pos2 " << pos2 << " distance " << distance << endl;
		cout << "vv (" << vel[0] << "," << vel[1] << "," << vel[2] << ")";
		cout << "gvv (" << gvv[0] << "," << gvv[1] << "," << gvv[2] << ") " << (distance - configuration.planetRadius) << " " << calcLatitude(pos) << " " << calcLongitude(pos) << endl;
		cout << "thrust2 " << thrust2 << " (" << controls[daeControlIndex[iphase-1].thrustX] << "," << controls[daeControlIndex[iphase-1].thrustY] << "," << controls[daeControlIndex[iphase-1].thrustZ] << ")" << endl;
		cout << "mu " << configuration.planetMu << " d^3 " << (pos2 * distance) << endl;
		cout << "drag " << dragFactor << " gravityFactor " << gravityFactor << " mass " << states[stateIndex[iphase-1].mass]<<endl;
		cout << "gravity " << sqrt(pos2) * gravityFactor;
		cout << " drag " << sqrt(dot(gvv,gvv,3)) * dragFactor;
		cout << " thrust " << sqrt(thrust2) * thrustFactor << endl;
		cout << "Accel(" << daeAccelX << "," << daeAccelY << "," << daeAccelZ<<")"<<endl;
		cout << "dm =" << -sqrt(thrust2) / (isp * constants.g0) * 1000<<endl;
*/
			};
}

void Launch::setupDae()
{
	int stages = configuration.stages.size();

	double timeToPitchover = 90;
	double altitudeBegin  = configuration.launchAltitude;
	double altitudeEnd;
//	double travelBegin = 0;
	double travelEnd;
	double velHoriBegin = initialState.velY;
	double velHoriEnd;
	double velVertBegin = 0;
	double velVertEnd;
	double alphaBegin = 0;
	double alphaEnd;
	double posXBegin = configuration.planetRadius + altitudeBegin;
	double posXEnd;
	double posYBegin = 0;
	double posYEnd;
	double velXBegin = 0;
	double velXEnd;
	double velYBegin = velHoriBegin;
	double velYEnd;
	double targetVelocity = sqrt(configuration.planetMu / (configuration.planetRadius + configuration.targetPeriapsis));
	for (int stage = 1;stage <= stages; stage++)
	{
		using PsoptInterface::GuessStruct;

		double timeBegin = timLower[stage-1];
		double timeEnd   = (stage == stages) ? configuration.timeToOrbit : timLower[stage];
		ControlStruct tS;
		tS.thrustX = psoptInterface.CreateNewControl(stage,-configuration.stages[stage-1].combinedThrust,configuration.stages[stage-1].combinedThrust,GuessStruct(timeEnd < timeToPitchover ? configuration.stages[stage-1].combinedThrust : 0));
		tS.thrustY = psoptInterface.CreateNewControl(stage,-configuration.stages[stage-1].combinedThrust,configuration.stages[stage-1].combinedThrust,GuessStruct(timeEnd < timeToPitchover ? 0 : configuration.stages[stage-1].combinedThrust));
		tS.thrustZ = psoptInterface.CreateNewControl(stage,-configuration.stages[stage-1].combinedThrust,configuration.stages[stage-1].combinedThrust,GuessStruct(0));
		daeControlIndex.push_back(tS);
		if (timeEnd < configuration.timeToOrbit) {
			altitudeEnd = tanh((timeEnd - (configuration.timeToOrbit / 2)) / configuration.timeToOrbit * 4) * configuration.targetPeriapsis / 2 + configuration.targetPeriapsis / 2;
			//cout << "temp " << ((4 * timeEnd - 2 * timeToOrbit)/timeToOrbit) <<endl;
			//cout << "temp " << (tanh((4 * timeEnd - 2 * timeToOrbit)/timeToOrbit)) <<endl;
			//cout << "temp " << (1 - pow(tanh((4 * timeEnd - 2 * timeToOrbit)/timeToOrbit),2)) <<endl;
			velVertEnd = 2 * configuration.targetPeriapsis / configuration.timeToOrbit * (1 - pow(tanh((4 * timeEnd - 2 * configuration.timeToOrbit)/configuration.timeToOrbit),2));
		} else {
			altitudeEnd = configuration.targetPeriapsis;
			velVertEnd = 0;
		}
		if (timeEnd < timeToPitchover) {
			travelEnd = 0;
			//cout << velHoriEnd << "vhe\n";
			velHoriEnd = velHoriBegin;
		} else if (timeEnd < configuration.timeToOrbit) {
			double temp = ((timeEnd - timeToPitchover) / (configuration.timeToOrbit - timeToPitchover));
			velHoriEnd = initialState.velY + temp * temp * (targetVelocity - initialState.velY);
			travelEnd += (timeEnd - timeBegin) * (velHoriEnd + velHoriBegin) / 2;
		} else {
			velHoriEnd = targetVelocity;
			travelEnd += targetVelocity * (timeEnd - timeBegin);
		}
		alphaEnd = travelEnd / (configuration.planetRadius + altitudeEnd);
		posXEnd = cos(alphaEnd) * (configuration.planetRadius + altitudeEnd);
		posYEnd = sin(alphaEnd) * (configuration.planetRadius + altitudeEnd);
		velXEnd = velVertEnd * cos(alphaEnd) - velHoriEnd * sin(alphaEnd);
		velYEnd = velVertEnd * sin(alphaEnd) + velHoriEnd * cos(alphaEnd);
		GuessStruct posXGuess = GuessStruct(configuration.planetRadius + altitudeBegin, configuration.planetRadius + altitudeEnd, alphaBegin+pi/2, alphaEnd+pi/2);
		GuessStruct posYGuess = GuessStruct(configuration.planetRadius + altitudeBegin, configuration.planetRadius + altitudeEnd, alphaBegin, alphaEnd);

		StateIndexStruct stSt;
		stSt.posX = psoptInterface.CreateNewState(stage, -configuration.planetSOI, configuration.planetSOI, posXGuess/*GuessStruct(posXBegin, posXEnd)*/);//, [this] () { return daeVel[0]; });
		stSt.posY = psoptInterface.CreateNewState(stage, -configuration.planetSOI, configuration.planetSOI, posYGuess/*GuessStruct(posYBegin, posYEnd)*/);//, [this] () { return daeVel[1]; });
		stSt.posZ = psoptInterface.CreateNewState(stage, -configuration.planetSOI, configuration.planetSOI, GuessStruct(initialState.posZ));//, [this] () { return daeVel[2]; });
		stSt.velX = psoptInterface.CreateNewState(stage, -configuration.shipMaxV, configuration.shipMaxV, GuessStruct(velXBegin, velXEnd));//, [this] () { return daeAccelX; });
		stSt.velY = psoptInterface.CreateNewState(stage, -configuration.shipMaxV, configuration.shipMaxV, GuessStruct(velYBegin, velYEnd));//, [this] () { return daeAccelY; });
		stSt.velZ = psoptInterface.CreateNewState(stage, -configuration.shipMaxV, configuration.shipMaxV, GuessStruct(0));//, [this] () { return daeAccelZ; });
		stSt.mass = psoptInterface.CreateNewState(stage, configuration.stages[stage-1].mass - configuration.stages[stage-1].fuel, configuration.stages[stage-1].mass, GuessStruct(configuration.stages[stage-1].mass, configuration.stages[stage-1].mass - configuration.stages[stage-1].fuel));
		stateIndex.push_back(stSt);
		std::cout << "stage " << stage << " tim(" << timeBegin << "," << timeEnd << ") " <<
				"posBegin(" << posXBegin << "," << posYBegin << ") " <<
				"posEnd(" << posXEnd << "," << posYEnd << ") " <<
				"velRelBegin(" << velVertBegin << "," << velHoriBegin << ") " <<
				"velRelEnd(" << velVertEnd << "," << velHoriEnd << ") " <<
				"velBegin(" << velXBegin << "," << velYBegin << ") " <<
				"velEnd(" << velXEnd << "," << velYEnd << ") " <<
				"alpha(" << alphaBegin << "," << alphaEnd << ")" << std::endl;

		daePathIndex.thrust2 = psoptInterface.CreateNewPath(stage, 0, configuration.stages[stage-1].combinedThrust * configuration.stages[stage-1].combinedThrust);
		daePathIndex.distance2 = psoptInterface.CreateNewPath(stage, configuration.planetRadius * configuration.planetRadius, configuration.planetSOI * configuration.planetSOI);
		daePathIndex.eccentricity2= psoptInterface.CreateNewPath(stage, 0 * 0, 0.995 * 0.995);

//		travelBegin = travelEnd;
		altitudeBegin = altitudeEnd;
		velHoriBegin = velHoriEnd;
		velVertBegin = velVertEnd;
		alphaBegin = alphaEnd;
		posXBegin = posXEnd;
		posYBegin = posYEnd;
		velXBegin = velXEnd;
		velYBegin = velYEnd;


		psoptInterface.SetDaeInitFunction(daeFunction);
	}

}
/**
 * altitude in m
 * ret in atm
 */
adouble Launch::calcPressure(adouble altitude)
{
	return configuration.planetP0 * exp(-altitude / configuration.planetScaleHeight);
}

void Launch::calcGroundVelocityVector(adouble* pos, adouble* vel, adouble* gvv)
{
  adouble pos_norm = sqrt(dot(pos, pos, 3));
  adouble latitude = calcLatitude(pos);
  adouble longitude = calcLongitude(pos);

  adouble factor = pos_norm * 2 * pi * cos(latitude) / configuration.planetRotationPeriod;

  gvv[0] = vel[0] + sin(longitude) * factor;
  gvv[1] = vel[1] - cos(longitude) * factor;
  gvv[2] = vel[2];

}

inline adouble Launch::calcLatitude(adouble* pos)
{
  return atan2(pos[2], sqrt(pos[0] * pos[0] + pos[1] * pos[1]));
}

inline adouble Launch::calcLongitude(adouble* pos)
{
  return atan2(pos[1], pos[0]);
}


/**
 * stage >= 0
 */
adouble Launch::calcIsp(adouble pressure, int stage)
{
	vector<Engine> * engines = & configuration.stages[stage].engines;
	adouble real_p = pressure;
	if (real_p > 1) {
		throw "can not handle pressures > 1 in calcIsp";
		real_p = 1;
	}
	adouble denominatorSum = 0;
	for (int i=0; i < (int)engines->size(); i++) {
		denominatorSum += engines->at(i).thrust / (engines->at(i).ispVac + real_p * (engines->at(i).isp1Atm - engines->at(i).ispVac));
	}
	adouble res = configuration.stages[stage].combinedThrust / denominatorSum;
	engines = NULL;
	return res;
}
adouble Launch::calcIsp(adouble pressure, adouble isp_0, adouble isp_vac)
{
  adouble real_p = pressure;
  if (real_p > 1) {
        real_p = 1;
  }
  return isp_0 * real_p + isp_vac * (1 - real_p);
}

void Launch::setupEvents()
{
	// initial state
	psoptInterface.CreateNewEvent(1, initialState.posX, initialState.posX, [this] (adouble* initial_states, adouble*, adouble*, adouble&, adouble&, adouble*)
	{ return initial_states[stateIndex[0].posX]; });
	psoptInterface.CreateNewEvent(1, initialState.posY, initialState.posY, [this] (adouble* initial_states, adouble*, adouble*, adouble&, adouble&, adouble*)
	{ return initial_states[stateIndex[0].posY]; });
	psoptInterface.CreateNewEvent(1, initialState.posZ, initialState.posZ, [this] (adouble* initial_states, adouble*, adouble*, adouble&, adouble&, adouble*)
	{ return initial_states[stateIndex[0].posZ]; });
	psoptInterface.CreateNewEvent(1, initialState.velX, initialState.velX, [this] (adouble* initial_states, adouble*, adouble*, adouble&, adouble&, adouble*)
	{ return initial_states[stateIndex[0].velX]; });
	psoptInterface.CreateNewEvent(1, initialState.velY, initialState.velY, [this] (adouble* initial_states, adouble*, adouble*, adouble&, adouble&, adouble*)
	{ return initial_states[stateIndex[0].velY]; });
	psoptInterface.CreateNewEvent(1, initialState.velZ, initialState.velZ, [this] (adouble* initial_states, adouble*, adouble*, adouble&, adouble&, adouble*)
	{ return initial_states[stateIndex[0].velZ]; });
	psoptInterface.CreateNewEvent(1, initialState.mass, initialState.mass, [this] (adouble* initial_states, adouble*, adouble*, adouble&, adouble&, adouble*)
	{ return initial_states[stateIndex[0].mass]; });

	// reaching periapsis
	psoptInterface.CreateNewEvent(configuration.stages.size(), configuration.planetRadius + configuration.targetPeriapsis, configuration.planetRadius + configuration.targetPeriapsis, [this] (adouble* , adouble* final_states, adouble*, adouble&, adouble&, adouble*) -> adouble
	{ return calcPeriapsis(final_states, stateIndex[configuration.stages.size()-1]); });

}

adouble Launch::calcPeriapsis(adouble* states, StateIndexStruct index)
{
	adouble p[3] = {states[index.posX], states[index.posY], states[index.posZ]};
	adouble v[3] = {states[index.velX], states[index.velY], states[index.velZ]};
	adouble h[3];
	cross(p, v, h);
	adouble ev[3];
	cross(v, h, ev);
	adouble posNorm = sqrt( dot(p, p, 3) );
	adouble velNorm2 = dot(v, v, 3);

	ev[0] = ev[0] / configuration.planetMu - p[0] / posNorm;
	ev[1] = ev[1] / configuration.planetMu - p[1] / posNorm;
	ev[2] = ev[2] / configuration.planetMu - p[2] / posNorm;
	adouble e_norm = sqrt(dot(ev, ev, 3));

	adouble a = 1 / (2 / posNorm - velNorm2 / configuration.planetMu);

	return a * (1 - e_norm);
}

void Launch::setupLinkages()
{
	psoptInterface.SetLinkagesInitFunction([this] (adouble* linkages, adouble* xad) -> void
			{
		adouble time_prev, time_next;
		adouble stat_prev[NR_OF_STATES], stat_next[NR_OF_STATES];
		adouble co_prev[NR_OF_CONTROLS], co_next[NR_OF_CONTROLS];

		int stages = configuration.stages.size();
		int index = 0;
		for(int iphase=1;iphase<stages;iphase++) {
			// get states
			time_prev = get_final_time(xad, iphase);
			time_next = get_initial_time(xad, iphase+1);
			get_final_states(stat_prev, xad, iphase);
			get_initial_states(stat_next, xad, iphase+1);
			get_final_controls(co_prev, xad, iphase);
			get_initial_controls(co_next, xad, iphase+1);

			// time
			linkages[index++] = time_prev - time_next;

			// position
			linkages[index++] = stat_prev[stateIndex[iphase-1].posX] - stat_next[stateIndex[iphase].posX];
			linkages[index++] = stat_prev[stateIndex[iphase-1].posY] - stat_next[stateIndex[iphase].posY];
			linkages[index++] = stat_prev[stateIndex[iphase-1].posZ] - stat_next[stateIndex[iphase].posZ];

			// velocity
			linkages[index++] = stat_prev[stateIndex[iphase-1].velX] - stat_next[stateIndex[iphase].velX];
			linkages[index++] = stat_prev[stateIndex[iphase-1].velY] - stat_next[stateIndex[iphase].velY];
			linkages[index++] = stat_prev[stateIndex[iphase-1].velZ] - stat_next[stateIndex[iphase].velZ];

			// mass
			adouble mass_difference = configuration.stages[iphase-1].mass - configuration.stages[iphase-1].fuel     - configuration.stages[iphase].mass;
			linkages[index++] = stat_prev[stateIndex[iphase-1].mass]-mass_difference-stat_next[stateIndex[iphase].mass];

			// pitch
			adouble prev_pos[3] = {stat_prev[stateIndex[iphase-1].posX], stat_prev[stateIndex[iphase-1].posY], stat_prev[stateIndex[iphase-1].posZ]};
			adouble next_pos[3] = {stat_next[stateIndex[iphase].posX],   stat_next[stateIndex[iphase].posY],   stat_next[stateIndex[iphase].posZ]};
			linkages[index++] = getDirectionPitchDestructively(prev_pos, co_prev) - getDirectionPitchDestructively(next_pos, co_next);

/*
			// orientation does not change drastically
			adouble prevControlSize = sqrt(dot(co_prev, co_prev,3));
			adouble nextControlSize = sqrt(dot(co_next, co_next,3));
			adouble prev_control;
			adouble next_control;
			adouble delta = 1;
			time_prev -= delta;
			time_next += delta;

			get_interpolated_control(&prev_control, 1, iphase, time_prev, xad);
			get_interpolated_control(&next_control, 1, iphase+1, time_next, xad);
			linkages[index++] = (prev_control -  co_prev[0]) / prevControlSize + (next_control - co_next[0]) / nextControlSize;
			get_interpolated_control(&prev_control, 2, iphase, time_prev, xad);
			get_interpolated_control(&next_control, 2, iphase+1, time_next, xad);
			linkages[index++] = (prev_control -  co_prev[1]) / prevControlSize + (next_control - co_next[1]) / nextControlSize;
			get_interpolated_control(&prev_control, 3, iphase, time_prev, xad);
			get_interpolated_control(&next_control, 3, iphase+1, time_next, xad);
			linkages[index++] = (prev_control -  co_prev[2]) / prevControlSize + (next_control - co_next[2]) / nextControlSize;
*/
/*
			get_control_derivative(&prev_control, 1, iphase, time_prev, xad);
			get_control_derivative(&next_control, 1, iphase+1, time_next, xad);
			linkages[index++] = prev_control/prevFactor - next_control/nextFactor;
			get_control_derivative(&prev_control, 2, iphase, time_prev, xad);
			get_control_derivative(&next_control, 2, iphase+1, time_next, xad);
			linkages[index++] = prev_control/prevFactor - next_control/nextFactor;
			get_control_derivative(&prev_control, 3, iphase, time_prev, xad);
			get_control_derivative(&next_control, 3, iphase+1, time_next, xad);
			linkages[index++] = prev_control/prevFactor - next_control/nextFactor;
*/
		}

/*
		// initial thrust
		get_initial_controls(co_prev, xad, 1);
		linkages[index++] = dot(co_prev, co_prev, 3) - configuration.stages[0].combinedThrust * configuration.stages[0].combinedThrust;
*/
			});
}
int Launch::calcLinkagesNumber()
{
	return (configuration.stages.size() - 1) * (NR_OF_STATES + 1 /* time */ + 1 /* pitch */ /* + 3 control derivatives */) /* + 1 initial thrust */;
}

/**
 * Calculate the angle between two vectors
 */
adouble Launch::getAngle(adouble* a, adouble* b)
{
	adouble dott = dot(a, b, 3);
	adouble sa = dot(a, a, 3);
	adouble sb = dot(b, b, 3);
	adouble den = sqrt(sa * sb);
	adouble c = dott / den;
	if (den > 0) {
		return acos(c);
	} else {
		return 0;
	}
}

/**
 * Arguments are changed
 */
adouble Launch::getDirectionPitchDestructively(adouble* pos, adouble* dir)
{
	adouble pos_norm = sqrt(dot(pos, pos, 3));
	adouble dir_norm = sqrt(dot(dir, dir, 3));
	pos[0] /= pos_norm;
	pos[1] /= pos_norm;
	pos[2] /= pos_norm;
	dir[0] /= dir_norm;
	dir[1] /= dir_norm;
	dir[2] /= dir_norm;

	return pi/2 - acos(dot(pos, dir, 3));
}

adouble Launch::getDirectionPitch(adouble* pos, adouble* dir)
{
	adouble pos_norm = sqrt(dot(pos, pos, 3));
	adouble dir_norm = sqrt(dot(dir, dir, 3));
	adouble pc[3] = {pos[0], pos[1], pos[2]};
	adouble dc[3] = {dir[0], dir[1], dir[2]};
	pc[0] /= pos_norm;
	pc[1] /= pos_norm;
	pc[2] /= pos_norm;
	dc[0] /= dir_norm;
	dc[1] /= dir_norm;
	dc[2] /= dir_norm;

	return pi/2 - acos(dot(pc, dc, 3));
}

void Launch::setupCosts()
{
	psoptInterface.SetEndpointCostFunction([this] (adouble* initial_states, adouble* final_states, adouble* parameters, adouble& t0, adouble& tf, adouble* xad, int iphase) -> adouble {
		if (iphase == (int)configuration.stages.size()) {
			return -final_states[stateIndex[iphase-1].mass];
		} else {
			return 0;
		}
	}

	);
}

void Launch::performPostprocessing(bool silent)
{
	DMatrix x, u, t;
	x = resultVectors.x;
	u = resultVectors.u;
	t = resultVectors.t;
	int stages = configuration.stages.size();
	long cols = t.GetNoCols();
	DMatrix pos = DMatrix(6, cols);
	DMatrix thr = DMatrix(4, cols);
	thr(colon(1,3),colon()) = u;
	thr(4, colon()) = Sqrt(sum(elemProduct(u,u)));
	pos(colon(1,3),colon()) = x(colon(1,3),colon());
	pos(4,colon()) = Sqrt(sum(elemProduct(x(colon(1,3),colon()), x(colon(1,3),colon()))));
	DMatrix vel = DMatrix(4, cols);
	vel(colon(1,3),colon()) = x(colon(4,6),colon());
	vel(4,colon()) =  Sqrt(sum(elemProduct(vel(colon(1,3),colon()), vel(colon(1,3),colon()))));
	DMatrix ecc = DMatrix(1, cols);
	DMatrix mass = x(colon(7,7),colon());
	DMatrix control = DMatrix(4, cols);
	DMatrix altitude = DMatrix(2,cols);
	altitude(1,colon()) = pos(4,colon()) - configuration.planetRadius;
	for(int i=1;i<=cols;i++) {
		adouble states[6];
		states[0] = pos(1,i);
		states[1] = pos(2,i);
		states[2] = pos(3,i);
		states[3] = vel(1,i);
		states[4] = vel(2,i);
		states[5] = vel(3,i);
		adouble p[3]; p[0] = states[0]; p[1] = states[1]; p[2] = states[2];
		adouble v[3]; v[0] = states[3]; v[1] = states[4]; v[2] = states[5];

		adouble ev[3];
		calcEccentricityVector(p, v, ev);
		adouble e_norm = sqrt(dot(ev, ev, 3));
		adouble pos_norm = sqrt(dot(p, p, 3));
		adouble vel_norm = sqrt(dot(v, v, 3));
		adouble a = (1 / (2 / (pos_norm) - (vel_norm) * (vel_norm) / configuration.planetMu));
		ecc(1, i) = e_norm.value();
	    pos(5, i) = (a * (1 - e_norm)).value();
	    pos(6, i) = (a * (1 + e_norm)).value();
	    int stage = static_cast<int>(ceil((double)i * stages / cols));
	    double maxthrust = configuration.stages[stage-1].combinedThrust;
	    control(1, i) = thr(4, i) / maxthrust * 100;
	    adouble p2[3] = { states[0], states[1], states[2] };
	    adouble thrv[3] = { thr(1,i), thr(2,i), thr(3, i) };
	    adouble gvv[3];
	    calcGroundVelocityVector(p2, v, gvv);
	    control(2, i) = getDirectionPitch(p2, thrv).getValue() / pi * 180;
	    control(3, i) = getDirectionPitchDestructively(p2, gvv).getValue() / pi * 180;
	    control(4, i) = control(2, i) - control(3, i);
	}
	altitude(2,colon()) = pos(6,colon()) - configuration.planetRadius;
	string title = configuration.name;
	bool show = ! silent;
	psoptInterface.plot_win(t,pos     ,title,"time(s)", "Distance from planet center (m)",  "x y z sum peri apo", "png", title + "-pos.png", show);
	psoptInterface.plot_win(t,vel     ,title,"time(s)", "Velocity (m/s)",  "x y z sum", "png", title + "-vel.png", show);
	psoptInterface.plot_win(t,thr     ,title,"time(s)", "Thrust (kN)",  "x y z sum", "png", title + "-thr.png", show);
	psoptInterface.plot_win(t,mass    ,title,"time(s)", "Mass (kg)", "mass", "png", title + "-mass.png", show);
	psoptInterface.plot_win(t,ecc     ,title,"time(s)", "Eccentricity", "ecce", "png", title + "-ecce.png", show);
	psoptInterface.plot_win(t,altitude,title,"time(s)", "Altitude (m) / Apoapsis (m)", "alt apo", "png", title + "-altitude.png", show);
	psoptInterface.plot_win(t,control ,title,"time(s)", "Control", "thrust(%) pitch(deg) ProGr(deg) AOA(deg)", "png", title + "-control.png", show);
	psoptInterface.plot_win(t,control ,title,"time(s)", "Control", "thrust(%) pitch(deg) ProGr(deg) AOA(deg)", "pdf", title + "-control.pdf", false);
}

void Launch::calcEccentricityVector(adouble* pos, adouble* vel, adouble* ev)
{
  adouble h[3];
  adouble pos_norm = sqrt(dot(pos, pos, 3));

  cross(pos, vel, h);
  cross(vel, h, ev);
  ev[0] = ev[0] / configuration.planetMu - pos[0] / pos_norm;
  ev[1] = ev[1] / configuration.planetMu - pos[1] / pos_norm;
  ev[2] = ev[2] / configuration.planetMu - pos[2] / pos_norm;

}

string Launch::getPitchThrust(vector<string> command)
{
	adouble tim = getArgAsDouble(1, command);
	adouble pos[3], vel[3], thr[3];

	string res = std::to_string(tim.value()) + " ";
	int nodes = length(resultVectors.t);

	linear_interpolation(&pos[0], tim, resultVectors.t, resultVectors.x(1, colon()), nodes);
	linear_interpolation(&pos[1], tim, resultVectors.t, resultVectors.x(2, colon()), nodes);
	linear_interpolation(&pos[2], tim, resultVectors.t, resultVectors.x(3, colon()), nodes);
	linear_interpolation(&vel[0], tim, resultVectors.t, resultVectors.x(4, colon()), nodes);
	linear_interpolation(&vel[1], tim, resultVectors.t, resultVectors.x(5, colon()), nodes);
	linear_interpolation(&vel[2], tim, resultVectors.t, resultVectors.x(6, colon()), nodes);
	linear_interpolation(&thr[0], tim, resultVectors.t, resultVectors.u(1, colon()), nodes);
	linear_interpolation(&thr[1], tim, resultVectors.t, resultVectors.u(2, colon()), nodes);
	linear_interpolation(&thr[2], tim, resultVectors.t, resultVectors.u(3, colon()), nodes);

	cout << "Thr (" << thr[0] << "," << thr[1] << "," << thr[2] << ")" << endl;
	int phase = psoptInterface.getPhaseOfTime(tim);

	adouble tr = sqrt(dot(thr, thr, 3));

	cout << tr << "= thrust, thrustmax = " << configuration.stages[getStageOfPhase(phase-1)].combinedThrust << endl;

	tr = tr / configuration.stages[getStageOfPhase(phase-1)].combinedThrust;

	res.append(std::to_string(tr.value()));
	adouble p = getDirectionPitch(pos, thr).getValue() / pi * 180;
	res.append(" ").append(std::to_string(p.value()));
	cout << "PT RES " << res << "\n";
	return res;
}

string Launch::getFinalPhaseTimes(vector<string> command)
{
	vector<adouble> endTimes = psoptInterface.getEndTimesOfPhases();
	string tim = "END_TIMES";
	for (int i = 0; i < (int)endTimes.size(); i++)
	{
		tim.append(" ").append(std::to_string(endTimes[i].value()));
	}
	return tim;
}

int Launch::getStageOfPhase(int phase)
{
	return phase-1;
}

void Launch::Test()
{
	if (false) {
		StateIndexStruct sis;
		stateIndex.push_back(sis);
		stateIndex[0].posX = 0;
		stateIndex[0].posY = 1;
		stateIndex[0].posZ = 2;
		stateIndex[0].velX = 3;
		stateIndex[0].velY = 4;
		stateIndex[0].velZ = 5;
		stateIndex[0].mass = 6;
		ControlStruct cs;
		daeControlIndex.push_back(cs);
		daeControlIndex[0].thrustX = 0;
		daeControlIndex[0].thrustY = 1;
		daeControlIndex[0].thrustZ = 2;
		daePathIndex.distance2 = 0;
		daePathIndex.thrust2 = 1;
		daePathIndex.eccentricity2 = 2;

		adouble states[7] = {630000, 1000, -1000, 300, 177, -1, 8000};
		adouble controls[3] = {200, 15, 1};
		adouble deri[7];
		adouble path[3];
		adouble filler[1];
		daeFunction(deri, path, states, controls, filler, *filler, filler, 1);
		printf("(%f,%f,%f)(%f,%f,%f) %f / %f : %f : %f",
				deri[0].value(), deri[1].value(), deri[2].value(),
				deri[3].value(), deri[4].value(), deri[5].value(), deri[6].value(),
				path[0].value(), path[1].value(), path[2].value());
	}



}

} // end namespace Aptor
