// *******************************************************************
// 60% transmittance surface
// *******************************************************************
void len60(string transmittanceFileName)
{
  ofstream os(transmittanceFileName.c_str());
  for(int ilambdaNM=300;ilambdaNM<=1200;ilambdaNM++){
    double lambda=ilambdaNM/1000.0;  //In micrometers
    std::cout<<"   0.00000     "<<lambda<<"    "<<"0.60000     0.00000"<<endl;
    os<<"   0.00000     "<<lambda<<"    "<<"0.60000     0.00000"<<endl;
  }
  return;
}
