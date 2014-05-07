/**
 * Socket.cpp - Socket interface to send commands to the optimizer
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

#include "Interface/Socket.h"

#include "Exception/ServerException.h"

#include <iostream>
#include <string>
#include <winsock2.h>

using std::cout;
using std::endl;
using std::string;

namespace Aptor {

/*
// s is a global static and not a class variable because including winsock2.h in the header leads to problems
static SOCKET listenSocket;
static sockaddr_in listenAddr;
static SOCKET writeSocket;
static sockaddr_in writeAddr;
static int writeSockaddrLen;
*/

struct Socket::WindowsPointers {
	SOCKET listenSocket = INVALID_SOCKET;
	sockaddr_in listenAddr;
	SOCKET writeSocket = INVALID_SOCKET;
	sockaddr_in writeAddr;
	int writeSockaddrLen = -1;
};


Socket::Socket(int _port):
	wp (new WindowsPointers()),
	port(_port)
{
}

Socket::~Socket() {}

void Socket::WaitForConnection()
{
	cout << "port " << port << endl;
    WSAData v;
	WSAStartup(MAKEWORD(2,2),&v);
	wp->listenSocket=socket(AF_INET,SOCK_STREAM,IPPROTO_TCP);
	if(wp->listenSocket==INVALID_SOCKET) {
		cout << "socket creation failed: " << WSAGetLastError() << endl;
		throw ServerException("Socket creation error");
	}
	wp->listenAddr.sin_family=AF_INET;
	wp->listenAddr.sin_addr.s_addr=inet_addr("127.0.0.1");
	wp->listenAddr.sin_port=htons(port);
	int res = bind(wp->listenSocket, (struct sockaddr*)&wp->listenAddr, sizeof(wp->listenAddr));
	if (res  == SOCKET_ERROR) {
		std::cout << "Error - when binding " << WSAGetLastError() << std::endl;
		closesocket(wp->listenSocket);
		WSACleanup();
		throw ServerException("Socket Binding Error");
	}
	res = listen(wp->listenSocket, 2);
	if (res == SOCKET_ERROR) {
		fprintf(stderr, "Listen failed: %d\n", WSAGetLastError());
		WSACleanup();
		throw ServerException("Socket Listening Error");
	}
	cout << "Socket Init done\n";

	cout << "Socket Waiting for connection\n";
	wp->writeSockaddrLen = sizeof(wp->writeAddr);
	wp->writeSocket = accept(wp->listenSocket, (struct sockaddr*)&wp->writeAddr, &wp->writeSockaddrLen);
	if (wp->writeSocket == INVALID_SOCKET) {
		fprintf(stderr, "Accept failed: %d\n", WSAGetLastError());
		closesocket(wp->listenSocket);
		WSACleanup();
		throw ServerException("Socket Accepting Error");
	}
	if (wp->writeAddr.sin_addr.S_un.S_un_b.s_b1 == 127 &&
			wp->writeAddr.sin_addr.S_un.S_un_b.s_b2 == 0 &&
			wp->writeAddr.sin_addr.S_un.S_un_b.s_b3 == 0 &&
			wp->writeAddr.sin_addr.S_un.S_un_b.s_b4 == 1) {
		cout << "Socket Connected to " << inet_ntoa(wp->writeAddr.sin_addr) << ":" << ntohs(wp->writeAddr.sin_port) << endl;
	} else {
		cout << "denied" << inet_ntoa(wp->writeAddr.sin_addr) << ":" << ntohs(wp->writeAddr.sin_port) << endl;
	}
}

string Socket::ReceiveMessage()
{
	string message = "";
	int transferred = BUFFLEN;
	while (transferred == BUFFLEN) {
		transferred = recv(wp->writeSocket,buffer,BUFFLEN, 0);
		if (transferred == SOCKET_ERROR) {
			cout << "Socket Error " << WSAGetLastError() << endl;
			throw ServerException("Socket Error");
		} else if (transferred == 0) {
			// socket closed
			return "SOCKET_CLOSED";
		} else {
			message.append(buffer, transferred);
		}
		cout << transferred << endl;
	}
	cout << message << endl;
	return message;
}

void Socket::WriteMessage(string message)
{
	cout << "OUT: " << message << endl;
	message.copy(buffer, message.size());
	send(wp->writeSocket, buffer, message.size(), 0);
}

void Socket::CloseConnection()
{
	closesocket(wp->listenSocket);
	WSACleanup();

}
void Socket::StopListening()
{

}


} // End namespace Aptor
