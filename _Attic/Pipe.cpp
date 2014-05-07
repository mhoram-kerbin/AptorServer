#include "Pipe.h"
#include <string>
#include <iostream>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace Aptor {

//#include "HandleObject.h"
using namespace std;

// ph is a global static and not a class variable because including windows.h in the header leads to problems
static HANDLE ph;


Pipe::Pipe(string name)
{
	string pipename = "\\\\.\\pipe\\" + name;
	cout << "creating pipe " << pipename << endl;
	ph = CreateNamedPipe(
			pipename.c_str(), // name
			PIPE_ACCESS_DUPLEX | // both directions
			FILE_FLAG_FIRST_PIPE_INSTANCE, // only one
			PIPE_TYPE_MESSAGE | PIPE_READMODE_MESSAGE | // message mode
			PIPE_WAIT | // blocking
			PIPE_REJECT_REMOTE_CLIENTS, // no remote clients
			1, // max instances
			4096, // out buffer size
			4096, // in buffer size
			0, // wait timeout
			NULL // Security Detail
	);
	if (ph == INVALID_HANDLE_VALUE)
	{
		printf("Could not create the pipe\n");
		exit(1);
	}
}

void Pipe::WaitForConnection()
{
	cout << "connecting ... ";
	ConnectNamedPipe(ph, NULL);
	cout << "DONE" << endl;
}

void Pipe::Write(string text)
{
	int size = text.size();
	if (size > BUFFER_SIZE) {
		throw "Write Command with too much data";
	}
	text.copy(buffer, size);
//	int written;
/*	if (WriteFile(
			ph, // handle
			buffer, // buffer to write
			size, // size of bytes to write
			&written, // number of written bytes
			NULL // not overlapped
			)) {
		cout << "write success:" << written << endl;
	} else {
		cout << "write error not handled" << endl;
	}
*/
}

string Pipe::ReadLine()
{
	return "";
}

} // End namespace Aptor
