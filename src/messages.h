// ************************************************************************* //
// Declaration of message services                                           //
//                                                                           //
// By Benjamin Fuks - 12.01.2022                                             //
// ************************************************************************* //


#ifndef MESSAGES_H
#define MESSAGES_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>     // Strings                        //
// ---------------------------------------------------- //

// ************************************************************************* //
// Definition of the message services class                                  //
// ************************************************************************* //
void InitMessages();
void debug  (const std::string&);
void info   (const std::string&);
void warning(const std::string&);
void error  (const std::string&);
void print(const std::string&, const std::string&);
void printerr(const std::string&);


#endif

