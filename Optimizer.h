/**
 * Optimizer.h - Base class of all optimizers
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

#include <ModernPsoptInterface.h>

#include <unordered_map>
#include <vector>
#include <string>
#include <functional>

using std::vector;
using std::string;
using std::exception;


namespace Aptor {


class Optimizer {
private:
	PsoptInterface::AlgoConfigStruct aCS;
protected:
	PsoptInterface::ModernPsoptInterface psoptInterface;
	std::unordered_map <string, std::function<string(vector<string>)> > commandMapper;
	vector<string> splitString(string command);
	double getArgAsDouble(int index, vector<string> command);
	int getArgAsInt(int index, vector<string> command);
	PsoptInterface::ResultVectors resultVectors;
	string getControlResults(vector<string> command);
	string getStateResults(vector<string> command);
	void initResults();
	void setIterations(vector<string> command);
	void setNodes(vector<string> command);
	void setNlpTolerance(vector<string> command);
	void setMeshRefinement(vector<string> command);
	void compute(string name, string output, int linkages);
public:
	Optimizer();
	string Command(string command);
};


} // end namespace Aptor
