#ifndef __MP_COMPLEX_H__
#define __MP_COMPLEX_H__

#include <string>
#include <iostream>
#include <mpc.h>

using namespace std;

class mpcomplex {
public:
    mpc_t mpc_val;
    mp_rnd_t   mpc_rnd;
    mp_prec_t  mpc_prec;


public:
    static mp_rnd_t   default_rnd;
    static mp_prec_t  default_prec;
    
    static mpcomplex PI();
    static mpcomplex PI(const mp_prec_t&, const mp_rnd_t&);
    static mpcomplex I();
    
    

public:
    // // Constructors && type conversion
    mpcomplex();
    mpcomplex( char* num , const mp_prec_t &p=default_prec, const mp_rnd_t &r=default_rnd );
    mpcomplex( const mpc_t& num , const mp_prec_t &p=default_prec, const mp_rnd_t &r=default_rnd );
    mpcomplex( const long double&, const long double& d=0 , const mp_prec_t &p=default_prec, const mp_rnd_t &r=default_rnd );
    mpcomplex( const mpcomplex& , const mp_prec_t &p=default_prec, const mp_rnd_t &r=default_rnd );
    mpcomplex( const mpfr_t& , const mp_prec_t &p=default_prec, const mp_rnd_t &r=default_rnd );
    mpcomplex( int i, const mpfr_t& imag,  const mp_prec_t &p=default_prec, const mp_rnd_t &r=default_rnd  );
    ~mpcomplex();

    mpcomplex operator=(const mpcomplex& o);

    string to_string(int precision=10);
    
    void init();
    void set_properties( mp_rnd_t r, mp_prec_t p );

    mpcomplex abs() const;
    mpcomplex neg() const;
    mpcomplex negate() const;
    mpcomplex& negate2();
    mpcomplex sqrt() const;
    void sqr3();
    mpcomplex& sqrt2();
    void sqrt3();
    mpcomplex square() const;
    mpcomplex& square2();
    mpcomplex Re() const; 
    mpcomplex Im() const;
    mpcomplex& timesI2();
    mpcomplex conju() const;
    mpcomplex conj2();
    mpcomplex log2();
    mpcomplex loga() const;
    mpcomplex exp2();
    mpcomplex expo() const;
    

    mpcomplex& operator+=( const mpcomplex& a);
    mpcomplex& operator-=( const mpcomplex& a);
    mpcomplex& operator*=( const mpcomplex& a);
    mpcomplex& operator/=( const mpcomplex& a);

    
    friend mpcomplex operator+(const mpcomplex& a, const mpcomplex& b);
    friend mpcomplex operator-(const mpcomplex& a, const mpcomplex& b);
    friend mpcomplex operator/(const mpcomplex& a, const mpcomplex& b);
    friend mpcomplex operator*(const mpcomplex& a, const mpcomplex& b);

    friend mpcomplex operator+(const mpcomplex& a, const long int& b);
    friend mpcomplex operator-(const mpcomplex& a, const long int& b);
    friend mpcomplex operator/(const mpcomplex& a, const long int& b);
    friend mpcomplex operator*(const mpcomplex& a, const long int& b);

    friend mpcomplex operator+(const long int& a, const mpcomplex& b);
    friend mpcomplex operator-(const long int& a, const mpcomplex& b);
    friend mpcomplex operator/(const long int& a, const mpcomplex& b);
    friend mpcomplex operator*(const long int& a, const mpcomplex& b);


    friend bool operator==(const mpcomplex& a, const mpcomplex& b);
    friend bool operator<(const mpcomplex& a, const mpcomplex& b);
    friend bool operator>(const mpcomplex& a, const mpcomplex& b);
    friend bool operator<=(const mpcomplex& a, const mpcomplex& b);
    friend bool operator>=(const mpcomplex& a, const mpcomplex& b);

};



#endif /* __MP_COMPLEX_H__ */
