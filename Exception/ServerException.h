/**
 * ServerException.h - Exception Class for the Aptor Server
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

#include <exception>
#include <string>
namespace Aptor {


class ServerException : public std::exception
{
private:
	std::string ex;
public:
	ServerException(std::string errStr) { ex = errStr; };
	const char* what() const noexcept { return ex.c_str(); };
};

} // end namespace Aptor
