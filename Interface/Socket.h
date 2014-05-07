/**
 * Socket.h - Socket interface to send commands to the optimizer
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

#include <string>
#include <memory>

#define BUFFLEN 1024

namespace Aptor {


class Socket {
private:
	struct WindowsPointers; // hide windows include
	std::unique_ptr<WindowsPointers> wp;

	char buffer[BUFFLEN];
	void ShutdownAfterSocketClosed();
	int port;
public:
	Socket(int port);
	~Socket();
	void WaitForConnection();
	std::string ReceiveMessage();
	void WriteMessage(std::string message);
	void CloseConnection();
	void StopListening();
};


} // end namespace Aptor
