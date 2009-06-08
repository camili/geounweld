
#include "mpcomplex.h"
#include <iostream>
#include <fstream>
#include <sstream> 
#include <vector>
#include <time.h>

using namespace std;

int vectors_size = 10000;
mpcomplex SQRT_2 = mpcomplex(2.0).sqrt2();
int output_precision = 200;

class GeounweldResults {
public:
    vector<mpcomplex> z;
    vector<mpcomplex> zcur;
    vector<mpcomplex> zpcur;
    vector<mpcomplex> b;
    vector<mpcomplex> a;
};

mpcomplex I(){
    return mpcomplex::I();
}
mpcomplex PI(){
    return mpcomplex::PI();
}

vector<mpcomplex> negateVector( vector<mpcomplex> points ) {
    vector<mpcomplex> other;
    for ( int i = 0 ; i < points.size(); i++) {
        other.push_back( points[i].negate() );
    }
    return other;
}

vector<mpcomplex> vectorFirstMap(vector<mpcomplex> z){
    vector<mpcomplex> result;
    result.push_back(z[0]);
    result.push_back(z[1]);
        for(int j = 2; j < z.size(); ++j)
        {                           
        result.push_back(((z[j]-z[1])/(z[j]-z[0])).sqrt2());
        }
    return result;
}
vector<mpcomplex> vectorZcur(vector<mpcomplex> z, int m){
    double dinf = 1e+20;
    vector<mpcomplex> result;
    result.push_back(mpcomplex(0.0,dinf));
    for(int j = 1; j < m; ++j) {
        mpcomplex zyd = z[0]+j*(z[1]-z[0])/m;      
        
        mpcomplex temp = ((zyd-z[1])/(zyd-z[0])).abs().sqrt2()*I();
        result.push_back(temp);
    }
    return result;
}

vector<mpcomplex> readFile( string file ) {
    vector<mpcomplex> numbers;

    string real, imag;
    ifstream file_op(file.c_str());

    while(!file_op.eof()) {
        real = "";
        imag = "";
        file_op >> real >> imag;
        if ( real.size() > 0 ) {
            string val = "(" + real + " " + imag + ")";
            char *stringCopy = new char[val.size()];
            strcpy(stringCopy,val.c_str());
            numbers.push_back( mpcomplex( stringCopy ) );
        }
    }         
    file_op.close();
    return numbers;
}

void writeArrayToFile(string file, vector<mpcomplex> array) {
    ofstream myfile;
    myfile.open(file.c_str());
    for ( int i = 0 ; i < array.size() ; i++ ) {
        myfile << array[i].to_string(output_precision) << endl;
    }
    myfile.close();
}     

void writeArrayToFile2(string file, vector<mpcomplex>& array1 , vector<mpcomplex>& array2 ) {
    ofstream myfile;
    myfile.open(file.c_str());
    for ( int i = 0 ; i < array1.size() ; i++ ) {
        myfile << array1[i].to_string(output_precision) << " " << array2[i].to_string(output_precision) << endl;
    }
    myfile.close();
}       

void writeArrayToFileBack(string file, vector<mpcomplex>& array) {
    ofstream myfile;
    myfile.open(file.c_str());
    for ( int i = 0 ; i < array.size() ; i++ ) {
        myfile << array[array.size()-1-i].to_string(output_precision) << endl;
    }
    myfile.close();
}

void printArray( vector<mpcomplex>& v, string label = "") {
    for(size_t i = 0; i < v.size(); ++i) {
        cout << label << i << "*" << v[i].to_string() << endl;
    }
}

inline mpcomplex calculateB( const mpcomplex& num ) {
    return num.Im()/num.Re();
}

vector<mpcomplex> calculateB( const vector<mpcomplex>& z, int n) {
    vector<mpcomplex> b;
    b.reserve(n);

    b.push_back(z[0].Im());
    b.push_back(z[1].Im());
    
    for(int i = 2; i < n; ++i) {
        b.push_back( calculateB(z[i]) );
    }
    return b;
}

inline mpcomplex calculateA( const mpcomplex& num, const mpcomplex& b ) {
    // return num.Re() + b*(num.Im());
    mpcomplex result = b;
    result*=num.Im();
    result+=num.Re();
    return result;
}

vector<mpcomplex> calculateA( const vector<mpcomplex>& z, const vector<mpcomplex>& b) {
    vector<mpcomplex> a;
    a.reserve(b.size());
    
    a.push_back(z[0].Re());
    a.push_back(z[1].Re());
    for(int i = 2; i < b.size(); ++i) {
        a.push_back( calculateA(z[i], b[i]) );
    }
    return a;
}

mpcomplex bodyOfLoop1( const mpcomplex& bj, const mpcomplex& aj, const mpcomplex& z) {
    mpcomplex yy = z.Im();
    mpcomplex y1 = yy;
    y1/=(aj-bj*yy);

    // mpcomplex y2 = (y1*y1+1.0).sqrt2();
    mpcomplex y2 = y1;
    y2.square2();
    y2+=1.0;
    y2.sqrt2();
    
    if (y1 < 0.0)
        return y2.negate2().timesI2();
    else
        return y2.timesI2();
}

inline mpcomplex bodyOfLoop2( const int& m , const int& k ) {
    // mpcomplex xdel = SQRT_2 - k*((SQRT_2-1.0)/m);
    mpcomplex xdel = SQRT_2;
    xdel-= k*((SQRT_2-1.0)/m);
    return (xdel.square2()-1.0).sqrt2().timesI2();
}

mpcomplex bodyOfLoop3( const mpcomplex& aj, const mpcomplex& zj, const mpcomplex& z ) {
    mpcomplex zz = z/(aj+zj*z);
    mpcomplex temp = zz.square2()-1.0;
    temp.sqrt2();
    if (temp.Im()*zz.Im() < 0.0)
        return temp.negate2();
    else
        return temp;
}


void recursiveZip( int m, int j , int n, 
    vector<mpcomplex>& z, vector<mpcomplex>& zcur, vector<mpcomplex>& zpcur, 
    vector<mpcomplex>& zNew, vector<mpcomplex>& zcurNew, vector<mpcomplex>& zpcurNew, 
    GeounweldResults &results ) {
        
    cout << "." ;
    cout.flush();
    mpcomplex bj = calculateB(z[j]);
    mpcomplex aj = calculateA(z[j],bj);
        
    results.b.push_back( bj );
    results.a.push_back( aj );
    
    mpcomplex zj = bj*I();

    for ( int i = 0 ; i < zcur.size() ; i++) {
        zcurNew.push_back(  bodyOfLoop1( bj,aj,zcur[i])  );
        zpcurNew.push_back( bodyOfLoop1( bj,aj,zpcur[i]) );
    }
    
    for ( int i = 0 ; i < m ; i++ ) {
        mpcomplex temp = bodyOfLoop2( m, i );
        zcurNew.push_back( temp );
        zpcurNew.push_back( temp.negate2() );
    }

    for ( int i = 0 ; i <= j ; i++ ) {
        zNew.push_back( z[i] );
    }
    for ( int i = j+1 ; i < z.size() ; i++ ) {
        zNew.push_back( bodyOfLoop3(aj, zj, z[i]) );
    }

    if ( j+1 < n ) {
        z.clear();
        zcur.clear();
        zpcur.clear();
        recursiveZip( m, j+1 , n, zNew, zcurNew, zpcurNew, z, zcur, zpcur, results);
    }
    else {
        results.zcur = zcurNew;
        results.z = zNew;
        results.zpcur = zpcurNew;
    }
}

void startRecursiveZip( vector<mpcomplex>& z, vector<mpcomplex>& zcur , const int& n, GeounweldResults &results) {
    vector<mpcomplex> zpcur = negateVector(zcur);
    
    vector<mpcomplex> zNew , zcurNew, zpcurNew;
    zNew.reserve( vectors_size );
    zcurNew.reserve(vectors_size);
    zpcurNew.reserve(vectors_size);
    
    // clock_t before = clock();

    recursiveZip( zcur.size() , 2, n , z , zcur, zpcur, zNew, zcurNew, zpcurNew,  results);

    // clock_t after = clock();
    
    cout << endl;
    // cout << "! " << (after - before)/(CLOCKS_PER_SEC/1000) << endl;
}

mpcomplex recalculatez(const mpcomplex& zcur, const mpcomplex& zc){
    mpcomplex zz = zcur/(mpcomplex(1.0)+ zc*zcur);
    return zz.square2();
}

vector<mpcomplex> recalculatezcur(const vector<mpcomplex>& zcur, const mpcomplex& zc){
    long double dinf = 1e+20;
    
    vector<mpcomplex> result;
    result.reserve(vectors_size);
    result.push_back(mpcomplex(-1*dinf, 0.0));
    
    for(int i = 1; i < zcur.size(); ++i) {
        mpcomplex temp = recalculatez(zcur[i], zc);
        result.push_back(temp);
    }
    
    result.push_back(mpcomplex(0.0, 0.0));
    return result;
}

vector<mpcomplex> recalculatez(const vector<mpcomplex>& z, const int& n, const mpcomplex& zc){
    vector<mpcomplex> result;
    result.reserve(vectors_size);
    
    for(size_t i = 0; i < n-1; ++i) {
        result.push_back(z[i]);
    }
    for(size_t i = n-1; i < z.size(); ++i) {
        mpcomplex temp = recalculatez(z[i], zc );
        result.push_back(temp);
    }
    return result;
}

vector<mpcomplex> finalzcur(const vector<mpcomplex>& zcur, const mpcomplex& znp){
    vector<mpcomplex> result;
    result.reserve(vectors_size);
    
    for(size_t i = 0; i < zcur.size(); ++i)
    {
        result.push_back((zcur[i]-znp)/(zcur[i]-znp.conju()));
    }
    return result;
}

mpcomplex logcomp(const mpcomplex& z){
    return z.negate().loga().Im()+PI();
}

vector<mpcomplex> logMap(const vector<mpcomplex>& zcur, const vector<mpcomplex>& zpcur){
    vector<mpcomplex> result;
    result.reserve(vectors_size);
    
    for(size_t i = 0; i < zcur.size(); ++i)
    {
        mpcomplex t = logcomp(zcur[i]) + logcomp(zpcur[i]).timesI2();
        result.push_back(t);
    }
    return result;
}     

void warningCheck( int i, const mpcomplex& zweld_i , const mpcomplex& xxo, const mpcomplex& yyo, const double& eps) {
    if((zweld_i.Re()-xxo).abs() < eps)  
    {
       int nxtest = 1;  
       int nxindx = i;    
       cout << " WARNING: data points for interior map are" << endl;
       cout << " too close in harmonic measure. Welding map" << endl;
       cout << " will have infinite slope. Index= " << nxindx << endl;
    }         
    if((zweld_i.Im()-yyo).abs() < eps) 
    {
       int nytest = 1;       
       int nyindx = i; 
       cout << " WARNING: data points for exterior map are" << endl;
       cout << " too close in harmonic measure. Welding map" << endl;
       cout << " will have slope = 0. Index= " << nyindx << endl;
    }
}

void iterativeWarningCheck( const vector<mpcomplex>& zweld ) {
    double eps = 1e-15;
    
    mpcomplex xxo = -1.0;
    mpcomplex yyo = -1.0;
    
    for ( int i = 0 ; i < zweld.size() ; ++i ) {
        // printf("%d/%d\n", i, zweld.size());
        warningCheck( i, zweld[i], xxo, yyo , eps);
        xxo = zweld[i].Re();
        yyo = zweld[i].Im();
    }
    
    // recursiveWarningCheck( 0 , zweld, mpcomplex(-1.0), mpcomplex(-1.0) );
}

void recursiveWarningCheck( int i , const vector<mpcomplex>& zweld, const mpcomplex& xxo, const mpcomplex& yyo ) {         
    // printf("%d/%d", i, zweld.size());
    double eps = 1e-15;
    warningCheck( i,zweld[i], xxo, yyo , eps);

    if ((i+1) < zweld.size()) recursiveWarningCheck( i+1, zweld, zweld[i].Re(), zweld[i].Im() );
}



int main(int argc, char* argv[]) {
    mpcomplex::default_prec = 100;
    mpcomplex::default_rnd = GMP_RNDU;
    
    string outputFolder = "results";
    string inputFile = "bdata.ga";
    
    // Data that can be entered by the user.
    int m = 5; // computed welding given at m*n points
    
    if ( argc > 1 ) {
        m = atoi(argv[1]);
    }
    if ( argc > 2 ) {
        mpcomplex::default_prec = atoi(argv[2]);
    }
    if ( argc > 3 ) {
        inputFile = argv[3];
    }
    if ( argc > 4 ) {
        outputFolder = argv[4];
    }

    system((string("mkdir -p ")+ outputFolder).c_str());
        
    // Reading and initializing data from file
    vector<mpcomplex> z = readFile(inputFile);
    
       
    mpcomplex zint = z.back();
    int n = z.size();
    
    // 
    vector<mpcomplex> zcur = vectorZcur(z, m);
    
    mpcomplex dinf1 = mpcomplex(1e-8);
    z.push_back(zint+dinf1);
    
     // Applying first map.
     vector<mpcomplex> map1z = vectorFirstMap(z);

     //Image of infinity after preliminary map.
     map1z.push_back(mpcomplex(1.0));
     map1z.push_back(mpcomplex(1.0 + dinf1));
     map1z.push_back(mpcomplex(1.0 - dinf1));
     
     vector<mpcomplex> b = calculateB(map1z, n);
     vector<mpcomplex> a = calculateA(map1z, b);
     
     for(size_t j = 2; j < n; ++j) {
         mpcomplex zj = b[j]*I();
     }
     
     GeounweldResults results;
     
     results.a.reserve(vectors_size);
     results.b.reserve(vectors_size);
     
     results.a.push_back( z[0].Re() );
     results.b.push_back( z[0].Im() );

     results.a.push_back( z[1].Re() );
     results.b.push_back( z[1].Im() );
     
     cout << "Number of points: " << n << endl;
     
     startRecursiveZip( map1z, zcur, n-1 , results);
     
     mpcomplex zc = mpcomplex(-1.0)/results.zcur[0];
     
     vector<mpcomplex> recalzcur = recalculatezcur(results.zcur, zc);
     vector<mpcomplex> recalzpcur = recalculatezcur(results.zpcur, zc);

     vector<mpcomplex> recalz = recalculatez(results.z, n, zc);
     
     mpcomplex znp = recalz[n-1]; 
     
     mpcomplex znpc = znp.conju();
     
     mpcomplex zfpz0 =  mpcomplex(1e+8)*((recalz[n]-znp)/(recalz[n]-znpc));   
    
     mpcomplex zpnp = recalz[n+1];   

     mpcomplex zpnpc = zpnp.conju();
     
     mpcomplex znp4 = (recalz[n+2]-zpnpc)/(recalz[n+2]-zpnp); 
     mpcomplex znp5 = (recalz[n+3]-zpnpc)/(recalz[n+3]-zpnp);  
     mpcomplex zainf = mpcomplex(2.0)*znp4*dinf1/(z[0]-z[1]);
     mpcomplex zbinf = (znp4+znp5)/mpcomplex(2.0)-(zainf*z[0]);    
     
     // Switch points if winding number is -1

     vector<mpcomplex> lastzcur, lastzpcur;
     if (znp.Im() < mpcomplex(0.0)) {
         lastzcur = finalzcur(recalzpcur, znp) ;
         lastzpcur = finalzcur(recalzcur, zpnpc);
     }
     else {
         lastzcur = finalzcur(recalzcur, znp) ;
         lastzpcur = finalzcur(recalzpcur, zpnpc);
     }

     vector<mpcomplex> zweld = logMap(lastzcur, lastzpcur);

     // Warn if points that are too close in harmonic measure. 

     recursiveWarningCheck( 0 , zweld, mpcomplex(-1.0), mpcomplex(-1.0) );
     // iterativeWarningCheck(zweld);
       
     if ( m != 1){
         mpcomplex argznn = (lastzcur.back()).loga().Im();  
         mpcomplex zd = ((argznn.negate()/m)*I()).expo();  
         mpcomplex argzpn = (lastzpcur.back()).loga().Im();  
         mpcomplex zpd = ((argzpn.negate()/m)*I()).expo();   
         for(size_t i = 1; i < m; ++i)
         {                                        
             lastzcur.push_back(lastzcur.back()*zd);   
             lastzpcur.push_back(lastzpcur.back()*zpd);    
             
             zweld.push_back(logcomp(lastzcur.back())+logcomp(lastzpcur.back())*I());
         }
         
     }             

     // After this some space is written in the files.
     //      write(4,*)'    '
     //      write(8,*)'    '
     //      x0=0.d0
     //      write(4,99)x0,x0
     //      write(8,99)x0,x0       
     
     if(znp.Im() > mpcomplex(0.0))
     {
         zweld.push_back(2*PI() + 2*PI()*I());
     }           
     else{
         zweld.push_back(mpcomplex(0.0, 0.0));
     }            


 // Need to reverse order if Im(znp)<0 so that welding starts at 0,0
 // and thus works with geoweld.
 
     if(znp.Im() < mpcomplex(0.0))
     {
         writeArrayToFileBack(outputFolder+"/welding-c.ga", zweld);       
     }  
     else{
         writeArrayToFile(outputFolder+"/welding-c.ga", zweld);        
     }            

     results.a.push_back(mpcomplex(0.0));     
     results.a.push_back(znp.Re());
     results.a.push_back(zpnp.Re());   
     results.a.push_back(z[0].Re());
     results.a.push_back(z[1].Re());
     results.a.push_back(zfpz0.Re());
     results.a.push_back(zainf.Re());
     results.a.push_back(zbinf.Re());     
     
     results.b.push_back(zc.Im());     
     results.b.push_back( znp.Im());
     results.b.push_back(zpnp.Im());   
     results.b.push_back(z[0].Im());
     results.b.push_back(z[1].Im());
     results.b.push_back(zfpz0.Im());
     results.b.push_back(zainf.Im());
     results.b.push_back(zbinf.Im());

     vector<mpcomplex> weldint;
     weldint.push_back(zint);
     weldint.push_back(z[0]);

     writeArrayToFile2(outputFolder+"/param-c.ga", results.a, results.b);
     writeArrayToFile(outputFolder+"/pverti-c.ga",lastzcur);
     writeArrayToFile(outputFolder+"/pverto-c.ga",lastzpcur);
     writeArrayToFile(outputFolder+"/weldinit-c.ga", weldint);
    return 0;
}