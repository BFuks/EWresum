// ************************************************************************* //
// Declaration of message services                                           //
//                                                                           //
// By Benjamin Fuks - 12.01.2022                                             //
// ************************************************************************* //


// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cstdlib>    // C Standard General Utilities   //
#include <fstream>    // Read/Write files               //
#include <iostream>   // In/Out streams                 //
// -------- Classes ----------------------------------- //
#include "messages.h" // Message services               //
// ---------------------------------------------------- //

// ************************************************************************* //
//  Initializing the log file                                                //
// ************************************************************************* //
void InitMessages()
{
   std::ofstream log("messages.log"); log.close();
}

// ************************************************************************* //
//  Printing information                                                     //
// ************************************************************************* //
void info(const std::string &str) { print(str,"INFO"); }
void warning(const std::string &str) { print(str,"WARNING");  }
void debug(const std::string &str) { print(str, "DEBUG");  }
void error(const std::string &str) { printerr(str); }

void print(const std::string &str, const std::string &tag)
{
   // log file
   std::ofstream log("messages.log", std::ios::app);
   log << " ** " << tag << ": " << str.c_str() << std::endl;
   log.close();

   // screen output
   if(tag.compare("WARNING")==0)
     std::cout << " ** \033[1;35mWARNING: " << str.c_str() <<"\033[0m"<< std::endl;
   else if(tag.compare("DEBUG")==0)
     std::cout << " ** \033[36mDEBUG  : " << str.c_str() <<"\033[0m"<< std::endl;
   else
     std::cout << " ** \033[33mINFO   : " << str.c_str() <<"\033[0m"<< std::endl;
}

void printerr(const std::string &str)
{
   std::cout << " ** \033[1;31mERROR  : " << str.c_str() <<"\033[0m"<< std::endl;
   std::cout << " ** \033[1;31mERROR  : Exiting... \033[0m"<< std::endl;
   std::ofstream log("messages.log", std::ios::app);
   log << " ** ERROR: " << str.c_str() << std::endl;
   log.close();
   exit(0);
}
