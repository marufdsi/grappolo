
 
#ifndef RNGSTREAM_H
#define RNGSTREAM_H
 
#include <string>

class RngStream
{
public:

RngStream (const char *name = "");


static bool SetPackageSeed (const unsigned comm_type seed[6]);


void ResetStartStream ();


void ResetStartSubstream ();


void ResetNextSubstream ();


void SetAntithetic (bool a);


void IncreasedPrecis (bool incp);


bool SetSeed (const unsigned comm_type seed[6]);


void AdvanceState (comm_type e, comm_type c);


void GetState (unsigned comm_type seed[6]) const;


void WriteState () const;


void WriteStateFull () const;


double RandU01 ();


int RandInt (int i, int j);



private:

double Cg[6], Bg[6], Ig[6];


bool anti, incPrec;


std::string name;


static double nextSeed[6];


double U01 ();


double U01d ();


};
 
#endif
 

