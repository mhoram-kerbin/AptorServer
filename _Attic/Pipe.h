#pragma once

#include <string>

namespace Aptor {

using std::string;

#define BUFFER_SIZE 4096

class Pipe
{
private:
	char buffer[BUFFER_SIZE];
public:
	Pipe(string name);
	void WaitForConnection();
	string ReadLine();
	void Write(string text);
};

} // End Namespace Aptor
