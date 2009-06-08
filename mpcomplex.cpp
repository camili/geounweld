#include "mpcomplex.h"

mp_rnd_t   mpcomplex::default_rnd = GMP_RNDD;
mp_prec_t  mpcomplex::default_prec = 50;



mpcomplex mpcomplex::PI() {
    return mpcomplex::PI(mpcomplex::default_prec,mpcomplex::default_rnd);
}

mpcomplex mpcomplex::I() {
    return mpcomplex(0.0, 1.0);
}

mpcomplex mpcomplex::PI(const mp_prec_t &p, const mp_rnd_t &r) {
    mpfr_t pi;
    mpfr_t one;
    mpfr_init2(one, p);
    mpfr_set_ui(one , 1, r);
    
    mpfr_init2(pi, p);
    mpfr_atan(pi, one , r );
    mpfr_mul_ui( pi , pi , 4 , r);
    return mpcomplex( pi );
}

mpcomplex::mpcomplex() { 
    set_properties( default_rnd, default_prec );
    init();
    mpc_set_ui(mpc_val,0, default_prec);
}

mpcomplex::mpcomplex( const mpc_t& num, const mp_prec_t &p, const mp_rnd_t &r ) {
    set_properties( r, p );
    init();
    mpc_set( mpc_val , num , mpc_prec );
}

mpcomplex::mpcomplex( const long double& real, const long double& imag, const mp_prec_t &p, const mp_rnd_t &r ) {
    set_properties( r, p );
    init();
    mpc_set_d_d( mpc_val , real, imag , mpc_prec );
}

mpcomplex::mpcomplex( const mpfr_t& real, const mp_prec_t &p, const mp_rnd_t &r ) {
    set_properties( r, p );
    init();
    mpc_set_fr(mpc_val, real, mpc_prec);
}

mpcomplex::mpcomplex( const mpcomplex &other, const mp_prec_t &p, const mp_rnd_t &r ) {
    set_properties( r, p );
    init();
    mpc_set( mpc_val , other.mpc_val , mpc_prec );
}

mpcomplex::mpcomplex( char* num, const mp_prec_t &p, const mp_rnd_t &r ) {
    set_properties( r, p );
    init();
    mpc_set_str(mpc_val, num , 10, mpc_rnd);
}

mpcomplex::mpcomplex( int i, const mpfr_t& imag, const mp_prec_t &p, const mp_rnd_t &r ) {
    set_properties( r, p );
    init();
    mpfr_t real;
    mpfr_init2(real, mpc_prec);
    mpfr_set_ui(real, i, mpc_rnd);
    mpc_set_fr_fr(mpc_val,  real, imag, mpc_rnd);
}

mpcomplex::~mpcomplex() { 
    mpc_clear(mpc_val); 
}

mpcomplex mpcomplex::operator=(const mpcomplex& o) {
    mpc_clear(mpc_val);
    init();
    mpc_set( mpc_val , o.mpc_val , mpc_prec );
    return *this;
}

void mpcomplex::init() {
    mpc_init3( mpc_val , mpc_prec, mpc_prec );
}

void mpcomplex::set_properties( mp_rnd_t r, mp_prec_t p ) {
    mpc_rnd = r;
    mpc_prec = p;
}

mpcomplex mpcomplex::abs() const{
    mpfr_t absoluteValue;
    mpfr_init2(absoluteValue, mpc_prec);
    mpc_abs(absoluteValue, mpc_val , default_rnd );
    return mpcomplex( absoluteValue );
}

mpcomplex mpcomplex::neg() const{
    mpc_t negativeValue;
    mpc_init3( negativeValue , mpc_prec, default_rnd );
    mpc_neg(negativeValue, mpc_val , default_rnd );
    return mpcomplex( negativeValue );
}

mpcomplex mpcomplex::negate() const{
    return mpcomplex(0.0)-*this;
}

mpcomplex& mpcomplex::negate2() {
    mpc_ui_sub( mpc_val , 0 , mpc_val , default_rnd);
    return *this;
}


mpcomplex mpcomplex::sqrt() const{
    mpc_t sqrtValue;
    mpc_init3(sqrtValue , mpc_prec, default_rnd );
    mpc_sqrt(sqrtValue, mpc_val , default_rnd );
    return mpcomplex( sqrtValue );
}

mpcomplex& mpcomplex::sqrt2() {
    mpc_sqrt(mpc_val, mpc_val , default_rnd );
    return *this;
}

void mpcomplex::sqrt3() {
    mpc_sqrt(mpc_val, mpc_val , default_rnd );
}

mpcomplex mpcomplex::square() const{
    return mpcomplex( *this ).square2();
}


mpcomplex& mpcomplex::square2() {
    mpc_mul(mpc_val, mpc_val, mpc_val, default_rnd);
    return *this;
}

void mpcomplex::sqr3() {
    mpc_sqr(mpc_val, mpc_val, default_rnd);
}


mpcomplex mpcomplex::Re() const {
    mpfr_t realPart;
    mpfr_init2(realPart, mpc_prec);
    mpc_real(realPart, mpc_val , default_rnd );
    return mpcomplex( realPart);
}

mpcomplex mpcomplex::Im() const{
    mpfr_t imagPart;
    mpfr_init2(imagPart, mpc_prec);
    mpc_imag(imagPart, mpc_val , default_rnd );
    return mpcomplex(imagPart);
}

mpcomplex& mpcomplex::timesI2() {
    return *this*=mpcomplex::I();
}

mpcomplex mpcomplex::conj2() {
    mpc_conj(mpc_val, mpc_val, default_rnd);
    return *this;
}
 
mpcomplex mpcomplex::conju() const{
    return mpcomplex( *this ).conj2();
}

mpcomplex mpcomplex::log2() {
    mpc_log(mpc_val, mpc_val, default_rnd);
    return *this;
}
 
mpcomplex mpcomplex::loga() const{
    return mpcomplex( *this ).log2();
}      

mpcomplex mpcomplex::exp2() {
    mpc_exp(mpc_val, mpc_val, default_rnd);
    return *this;
}
 
mpcomplex mpcomplex::expo() const{
    return mpcomplex( *this ).exp2();
}

mpcomplex& mpcomplex::operator+=( const mpcomplex& a) {
    mpc_add(mpc_val, mpc_val, a.mpc_val, default_rnd);
    return *this;
}

mpcomplex& mpcomplex::operator-=( const mpcomplex& a) {
    mpc_sub(mpc_val, mpc_val, a.mpc_val, default_rnd);
    return *this;
}

mpcomplex& mpcomplex::operator*=( const mpcomplex& a) {
    mpc_mul(mpc_val, mpc_val, a.mpc_val, default_rnd);
    return *this;
}

mpcomplex& mpcomplex::operator/=( const mpcomplex& a) {
    mpc_div(mpc_val, mpc_val, a.mpc_val, default_rnd);
    return *this;
}


/*
 * These functions operate on complex and complex
 */
mpcomplex operator+(const mpcomplex& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_add(value, a.mpc_val, b.mpc_val, a.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator-(const mpcomplex& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_sub(value, a.mpc_val, b.mpc_val, a.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator*(const mpcomplex& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_mul(value, a.mpc_val, b.mpc_val, a.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator/(const mpcomplex& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_div(value, a.mpc_val, b.mpc_val, a.default_rnd);
    return mpcomplex(value);
}


/*
 * These functions operate on complex and long integers
 */
mpcomplex operator+(const mpcomplex& a, const long int& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_add_ui(value, a.mpc_val, b, a.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator-(const mpcomplex& a, const long int& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_sub_ui(value, a.mpc_val, b, a.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator*(const mpcomplex& a, const long int& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_mul_ui(value, a.mpc_val, b, a.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator/(const mpcomplex& a, const long int& b) {
    mpc_t value;
    mpc_init3( value , a.mpc_prec, a.mpc_prec );
    mpc_div_ui(value, a.mpc_val, b, a.default_rnd);
    return mpcomplex(value);
}


mpcomplex operator+(const long int& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , b.mpc_prec, b.mpc_prec );
    mpc_add_ui(value, b.mpc_val, a, b.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator-(const long int& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , b.mpc_prec, b.mpc_prec );
    mpc_ui_sub(value, a , b.mpc_val, b.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator*(const long int& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , b.mpc_prec, b.mpc_prec );
    mpc_mul_ui(value, b.mpc_val, a, b.default_rnd);
    return mpcomplex(value);
}

mpcomplex operator/(const long int& a, const mpcomplex& b) {
    mpc_t value;
    mpc_init3( value , b.mpc_prec, b.mpc_prec );
    mpc_ui_div(value, a, b.mpc_val, b.mpc_rnd);
    return mpcomplex(value);
}

/*
 * These functions operate on complex and long doubles
 */


bool operator==(const mpcomplex& a, const mpcomplex& b) {
    return mpc_cmp(a.mpc_val,b.mpc_val) == 0;
}

bool operator<(const mpcomplex& a, const mpcomplex& b) {
    int result = mpc_cmp(a.mpc_val,b.mpc_val);
    return MPC_INEX_RE(result)<0;
}

bool operator>(const mpcomplex& a, const mpcomplex& b) {
    int result = mpc_cmp(a.mpc_val,b.mpc_val);
    return MPC_INEX_RE(result)>0;
}

bool operator<=(const mpcomplex& a, const mpcomplex& b) {
    int result = mpc_cmp(a.mpc_val,b.mpc_val);
    return MPC_INEX_RE(result)<=0;
}

bool operator>=(const mpcomplex& a, const mpcomplex& b) {
    int result = mpc_cmp(a.mpc_val,b.mpc_val);
    return MPC_INEX_RE(result)>=0;
}

string mpcomplex::to_string(int precision) {
    int base = 10;
    char * stringRepresentation = mpc_get_str( base , precision , mpc_val , mpc_rnd );
    string returnValue = stringRepresentation;
    mpc_free_str( stringRepresentation );
    return returnValue;
}
