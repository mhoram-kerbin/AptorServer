/**
 * Launch.h - A Launch Trajectory Optimizer for KSP based on PSOPT
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

#pragma once

#include "Optimizer.h"

#include <functional>

namespace Aptor {

using std::string;

class Launch : public Optimizer
{
private:
	struct Engine {
		double thrust;  // thrust in kN
		double isp1Atm; // ISP at Kerbin Seavelevel
		double ispVac;  // ISP in vacuum
		Engine(double _thrust, double _isp1Atm, double _ispVac)
		{
			thrust = _thrust;
			isp1Atm = _isp1Atm;
			ispVac = _ispVac;
		}
	};
	struct StageParameters {
		double mass;    // 1. total mass of ship at begin of the stage in kg
		double fuel;    // 2. mass of fuel in this stage in kg
		vector<Engine> engines; // 4. Active engines in this stage
		double combinedThrust = 0; // thrust based on the engines
		double drag;    // 6. drag coefficient (dimensionless)
	};
	struct {
		string name                 = "Launch Optimization";
		double planetMass           = 0; // in kg
		double planetRadius         = 0; // in m
		double planetMu             = 0; // in m^3 s^-2
		double planetScaleHeight    = 0; // in m
		double planetP0             = 0; // in atm
		double planetRotationPeriod = 0; // in s
		double planetSOI            = 0; // in m
		double launchLatitude       = 0; // in rad
		double launchLongitude      = 0; // in rad
		double launchAltitude       = 0; // in m
		double shipMaxV             = 0; // in m s^-1
		double targetPeriapsis      = 0; // in m
		double timeToOrbit          = 360; // in s
		vector<StageParameters> stages;
	} configuration;
	const struct {
		double conversionFactor = 1.2230948554874; // conversion factor for pressure and density in kg m^-3 atm^-1
		double gravitationalConstant = 6.674E-11; // gravitational constant in N m kg^-2
		double g0 = 9.82; // conversion constant for fuel consumption of engines in m s^-2
	} constants;
	struct StateIndexStruct {int posX; int posY; int posZ; int velX; int velY; int velZ; int mass;};
	vector<StateIndexStruct> stateIndex;
	adouble daeAccelX, daeAccelY, daeAccelZ;
	struct ControlStruct {int thrustX; int thrustY; int thrustZ;};
	vector<ControlStruct> daeControlIndex;
	struct {int thrust2; int distance2; int eccentricity2;} daePathIndex;
	vector<double> timLower;
	vector<double> timUpper;
	struct {
		double posX;
		double posY;
		double posZ;
		double velX;
		double velY;
		double velZ;
		double mass;
	} initialState;

	std::function<void (adouble* derivatives, adouble* path, adouble* states,
					adouble* controls, adouble* parameters, adouble& time,
					adouble* xad, int iphase)> daeFunction;
	void setupDaeFunction();

	void setName(vector<string> command);
	void setPlanetMass(vector<string> command);
	void setPlanetRadius(vector<string> command);
	void setPlanetScaleHeight(vector<string> command);
	void setPlanetP0(vector<string> command);
	void setPlanetRotationPeriod(vector<string> command);
	void setPlanetSOI(vector<string> command);
	void setLaunchLatitude(vector<string> command);
	void setLaunchLongitude(vector<string> command);
	void setLaunchAltitude(vector<string> command);
	void setStages(vector<string> command);
	void setShipMaxV(vector<string> command);
	void addStage(vector<string> command);
	void addEngine(vector<string> command);
	void setTargetPeriapsis(vector<string> command);
	void compute(vector<string> command);

	adouble calcPressure(adouble altitude);
	void calcGroundVelocityVector(adouble* pos, adouble* vel, adouble* gvv);
	inline adouble calcLatitude(adouble* pos);
	inline adouble calcLongitude(adouble* pos);

	adouble calcIsp(adouble pressure, int stage);
	adouble calcIsp(adouble pressure, adouble isp_0, adouble isp_vac);
	adouble calcPeriapsis(adouble* states, StateIndexStruct index);
	int calcLinkagesNumber();
	void setupInitialState();
	void setupTime();
	void setupDae();
	void setupEvents();
	void setupLinkages();
	void setupCosts();
	void performPostprocessing(bool silent = true);
	void calcEccentricityVector(adouble* pos, adouble* vel, adouble* ev);
	int getStageOfPhase(int phase);
	adouble getAngle(adouble* a, adouble* b);

	adouble getDirectionPitchDestructively(adouble* pos, adouble* dir);
	adouble getDirectionPitch(adouble* pos, adouble* dir);

protected:
	string getPitchThrust(vector<string> command);
	string getFinalPhaseTimes(vector<string> command);
public:
	Launch();
	void Test();
};

} // end namespace Aptor
