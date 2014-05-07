/**
 * AptorServer.cxx - Main Entry for Aptor Server
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

#include "Control.h"
#include "Exception/ServerException.h"
#include "Optimizer/Launch.h"
#include "Interface/Socket.h"

#include <cstdio>
#include <sstream>
#include <string>
#include <fstream>

namespace Aptor {

using namespace std;

void readFromFile(string file);
void readFromPipe(string name);
void readFromSocket(string name);

void readFromFile(string file)
{
	std::ifstream conf(file);
	std::string line;

	Launch l;
	while (std::getline(conf, line))
	{
		try {
			cout << l.Command(line) << endl;
		} catch (const ServerException& e) {
			cout << e.what() << endl;
		}
	}
	conf.close();
}

void readFromPipe(string name)
{
/*
	Pipe p = Pipe(name);
	p.WaitForConnection();
	p.Write("Aptor Server V" + to_string(APTOR_VERSION) + " accepting connections");
	string res = p.ReadLine();
	cout << "Got Message: \"" + res + "\"" << endl;
*/
}

void readFromSocket(string name)
{
	Launch l;
	int port = stoi(name);
	Socket s(port);
	bool waitForConnection = true;
	while (waitForConnection) {
		s.WaitForConnection();
		bool receiveMessages = true;
		while (receiveMessages) {
			string command = s.ReceiveMessage();
			if (command == "SOCKET_CLOSED") {
				receiveMessages = false;
			} else if (command == "SHUTDOWN") {
				receiveMessages = false;
				waitForConnection = false;
			} else {
				string res;
				try {
					res = l.Command(command);
				} catch (const ServerException & e) {
					res = "FAILURE: ";
					res.append(e.what());
				}
				s.WriteMessage(res);
			}
		}
		s.CloseConnection();
	}
//	s.SendHello();
//	p.WaitForConnection();
//	p.Write("Aptor Server V" + to_string(APTOR_VERSION) + " accepting connections");
//	string res = p.ReadLine();
//	cout << "Got Message: \"" + res + "\"" << endl;

}

} // end namespace Aptor

using namespace Aptor;

int main(int argc, char* argv[])
{
	vector<string> _argv;
	for (int i = 0; i < argc; i++) {
		_argv.push_back(argv[i]);
	}
	Control c(argc, _argv);
	c.Run();

	if(argc > 1) {
		string c = argv[1];
		if (c == "file") {
			if (argc > 2) {
				readFromFile(argv[2]);
			} else {
				readFromFile("config.h");
			}
		} else if (c == "pipe") {
			if (argc > 2) {
				readFromPipe(argv[2]);
			} else {
				readFromPipe("AptorServer");
			}
		} else if (c == "socket") {
			if (argc > 2) {
				readFromSocket(argv[2]);
			} else {
				readFromSocket("AptorServer");
			}
		}
	} else {
		readFromFile("config.h");
	}

	Launch l();

//	l.Test();

	return 1;
}

