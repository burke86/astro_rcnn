// *******************************************************************
// 40% refklectance mirror
// *******************************************************************
void al40(string reflectanceFileName)
{
  ofstream os(reflectanceFileName.c_str());
  for(int ilambdaNM=300;ilambdaNM<=1200;ilambdaNM++){
    double lambda=ilambdaNM/1000.0;  //In micrometers
    std::cout<<"   0.00000     "<<lambda<<"    "<<"0.40000     0.00000"<<endl;
    os<<"   0.00000     "<<lambda<<"    "<<"0.40000     0.00000"<<endl;
  }
  return;
}
