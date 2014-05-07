/**
 * Optimizer.cpp - Base class of all optimizers
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

#include "Optimizer.h"

#include "Exception/ServerException.h"

#include <psopt.h>
#include <iostream>
#include <sstream>

namespace Aptor {

Optimizer::Optimizer()
{
	commandMapper["ITERATIONS"]             = [this] (vector<string> c) -> string { setIterations(c);     return "OK"; };
	commandMapper["SET_NODES"]              = [this] (vector<string> c) -> string { setNodes(c);          return "OK"; };
	commandMapper["NLP_TOLERANCE"]          = [this] (vector<string> c) -> string { setNlpTolerance(c);   return "OK"; };
	commandMapper["MESH_REFINEMENT"]        = [this] (vector<string> c) -> string { setMeshRefinement(c); return "OK"; };
	commandMapper["GET_CONTROLS"] = [this] (vector<string> c) -> string { return getControlResults(c); };
}

void Optimizer::setIterations(vector<string> command)
{
	aCS.nlpIterMax = getArgAsInt(1, command);
}
void Optimizer::setNodes(vector<string> command)
{
	string res = "[" + std::to_string(getArgAsInt(1,command));

	for (int i = 2;i< (int)command.size();i++) {
		res = res + "," + std::to_string(getArgAsInt(i,command));
	}
	res += "]";
	aCS.nodes = res;
}
void Optimizer::setNlpTolerance(vector<string> command)
{
	aCS.nlpTolerance = getArgAsDouble(1, command);
}
void Optimizer::setMeshRefinement(vector<string> command)
{
	aCS.mesh_refinement = command[1];
}
void Optimizer::compute(string name, string output, int linkages)
{
	psoptInterface.Psopt(name, output, linkages, aCS);
}
string Optimizer::Command(string command)
{
	if (command.size() == 0) {
		return "Empty line ignored";
	}
	if (command[0] == '#') {
		return "Comment ignored";
	}
	vector<string> tokens = splitString(command);
	string res;
	try {
		res = commandMapper.at(tokens[0])(tokens);
	} catch (const std::out_of_range& e) {
		throw ServerException("Command unknown :\"" + command + "\"");
	}
	return res;
}

vector<string> Optimizer::splitString(string command)
{
	string temp;
	std::stringstream stream(command);
	vector<string> tokens;
	while (stream >> temp) {
		tokens.push_back(temp);
	}
	return tokens;
}

/**
 * index > 0
 */
double Optimizer::getArgAsDouble(int index, vector<string> command)
{
	double ret;
	if (command.size() > 1) {
		try {
			ret = stod(command[index]);
		} catch (const std::invalid_argument& e) {
			throw ServerException("Argument to \"" + command[0] + "\" is not a number");
		}catch (const std::out_of_range& e) {
			throw ServerException("Argument to \"" + command[0] + "\" is out of range");
		}
	} else {
		throw ServerException("Argument missing");
	}
	return ret;
}

int Optimizer::getArgAsInt(int index, vector<string> command)
{
	int ret;
	if (command.size() > 1) {
		try {
			ret = stoi(command[index]);
		} catch (const std::invalid_argument& e) {
			throw ServerException("Argument to \"" + command[0] + "\" is not a number");
		}catch (const std::out_of_range& e) {
			throw ServerException("Argument to \"" + command[0] + "\" is out of range");
		}
	} else {
		throw ServerException("Argument missing for Command");
	}
	return ret;
}

void Optimizer::initResults()
{
	resultVectors = psoptInterface.getResults();
}

string Optimizer::getControlResults(vector<string> command)
{
	double startTime = getArgAsDouble(1, command);
	double endTime = startTime;
	if (command.size() > 2) {
		endTime = getArgAsDouble(2, command);
	}
	double stepTime = 1;
	if (command.size() > 3) {
		stepTime = getArgAsDouble(3, command);
	}
	//resultVectors.t.Print("T");
	vector<string> res;
	adouble value = 0;
	int nodes = resultVectors.t.GetNoCols();
	int controls = resultVectors.u.GetNoRows();
	adouble time;
	for (time = startTime; time <= endTime; time += stepTime) {
		res.push_back(std::to_string(time.value()));
	}


	for (int i = 1;i <= controls; i++) {
		//cout << "row " << i << endl;
		DMatrix tempU = resultVectors.u(i, colon());
		int counter = 0;
		for (time = startTime; time <= endTime; time += stepTime) {
			//cout << time << endl;
			linear_interpolation(&value, time, resultVectors.t, tempU, nodes);
			res[counter++] += " " + std::to_string(value.value());
		}
	}
	string ret = "";
	for (int i = 0; i < (int)res.size(); i++) {
		ret += res[i] + "\n";
	}

	return ret;
}
string Optimizer::getStateResults(vector<string> command)
{
	double startTime = getArgAsDouble(1, command);
	double endTime = startTime;
	if (command.size() > 2) {
		endTime = getArgAsDouble(2, command);
	}
	double stepTime = 1;
	if (command.size() > 3) {
		stepTime = getArgAsDouble(3, command);
	}
	//resultVectors.t.Print("T");
	vector<string> res;
	adouble value = 0;
	int nodes = resultVectors.t.GetNoCols();
	int states = resultVectors.x.GetNoRows();
	adouble time;
	for (time = startTime; time <= endTime; time += stepTime) {
		res.push_back(std::to_string(time.value()));
	}


	for (int i = 1;i <= states; i++) {
		//cout << "row " << i << endl;
		DMatrix tempX = resultVectors.x(i, colon());
		int counter = 0;
		for (time = startTime; time <= endTime; time += stepTime) {
			//cout << time << endl;
			linear_interpolation(&value, time, resultVectors.t, tempX, nodes);
			res[counter++] += " " + std::to_string(value.value());
		}
	}
	string ret = "";
	for (int i = 0; i < (int)res.size(); i++) {
		ret += res[i] + "\n";
	}

	return ret;
}


} // namespace Aptor
