// *******************************************************************
// From: http://refractiveindex.info/?group=GLASSES&material=BK7
// *******************************************************************
void BK7(string dispersionFileName)
{
  ofstream os(dispersionFileName.c_str());
  for(int ilambdaNM=300;ilambdaNM<=1530;ilambdaNM++){
    double lambda=ilambdaNM/1000.0;  //In micrometers
    double lambda2=lambda*lambda;
    double n2=(  (1.03961212*lambda2/(lambda2-0.00600069867)) + 
		 (0.231792344*lambda2/(lambda2-0.0200179144)) +
		 (1.01046945*lambda2/(lambda2-103.560653))    +
		 1.0);
    double n=sqrt(n2);
    std::cout<<lambda<<" "<<n<<std::endl;
    os<<lambda<<" "<<n<<std::endl;
  }
  return;
}
