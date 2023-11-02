// Full6DOF.cpp : Defines the entry point for the application.
//

#include "Full6DOF.h"

#include "lib/tinyxml/tinyxml2.h"

#include <iostream>
#include <cstring>
#include <string>

void display_help()
{

}

void run(char** argv)
{
	tinyxml2::XMLDocument simDocument;
	auto err = simDocument.LoadFile(argv[2]);
	if (err != tinyxml2::XML_SUCCESS) 
	{
		//Could not load file. Handle appropriately.
		throw std::invalid_argument("Could not load file.");
	}
}

void test(char** argv)
{

}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		std::cerr << "please provide a command." << std::endl;
		return -1;
	}
	if (argv[1][0] != '-' || strlen(argv[1]) != 2)
	{
		std::cerr << "please input a valid command." << std::endl;
	}
	switch (argv[1][1])
	{
	case 'h':
		display_help();
		break;
	case 'r':
		run(argv);
		break;
	case 't':
		test(argv);
		break;
	default:
		std::cerr << "invalid command." << std::endl;
		display_help();
		break;
	}

	return 0;
}
