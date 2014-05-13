// Ask libo for details.

#ifndef CAMERA_H_
#define CAMERA_H_

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

namespace idl {
namespace geometry {

template<typename Derived> class Traits;

template<typename T, int DataSize> class OwnedStorage;
template<typename T, int DataSize> class SharedStorage;

template<template<typename, int> class _Storage>
class CheckOwnedStorage;
template<>
class CheckOwnedStorage<OwnedStorage> {};

/** 
 * \brief A base class whose inheritants will have their own storage of data.
 */
template<typename T, int DataSize>
class OwnedStorage
{
private:
    T data_[DataSize];
public:
    OwnedStorage() {}
    explicit OwnedStorage(const T* data) { memcpy(data_, data, sizeof(T) * DataSize); }
    OwnedStorage(const OwnedStorage<T, DataSize>& other) { memcpy(data_, other.data(), sizeof(T) * DataSize); }
    OwnedStorage(const SharedStorage<T, DataSize>& other) { memcpy(data_, other.data(), sizeof(T) * DataSize); }
    
    inline T* data() { return data_; }
    inline const T* data() const { return data_; }
    inline int dataSize() const { return DataSize; }
    
    //OwnedStorage& operator =(const OwnedStorage& other) { memcpy(data_, other.data_, sizeof(T) * DataSize); return *this; }

protected:
    ~OwnedStorage() {}
};

/** 
 * \brief A base class whose inheritants will share the data buffer provided.
 */
template<typename T, int DataSize>
class SharedStorage 
{
private:
    T* data_;
    SharedStorage() {}
    SharedStorage(const SharedStorage& other) {}
public:
    explicit SharedStorage(T* data) : data_(data) {}
    
    inline T* data() { return data_; }
    inline const T* data() const { return data_; }
    inline int dataSize() const { return DataSize; }

    SharedStorage& operator =(const SharedStorage& other) { memcpy(data_, other.data_, sizeof(T) * DataSize); return *this; }

protected:
    ~SharedStorage() {}
};

template<typename Distortion,
    template<typename, int> class _Storage = OwnedStorage>
class Camera : public _Storage<typename Traits<Camera<Distortion, _Storage> >::Scalar, Traits<Camera<Distortion, _Storage> >::kDataSize>
{
    typedef Camera<Distortion, _Storage> Self;
    typedef typename Traits<Self>::Storage Storage;
    typedef typename Traits<Self>::Scalar T;
public:
    Camera();
    Camera(T* data);
    Camera(const T* data);
    Camera(const Self& other) : Storage(other) {}

public:
    void project(const T incident[3], T pixel[2], T jacobRay[3], T* jacobParam) const;
    const T& fx() const { return Storage::data()[0]; }
    const T& fy() const { return Storage::data()[1]; }
    const T& px() const { return Storage::data()[2]; }
    const T& py() const { return Storage::data()[3]; }
    const T& skew() const { return Storage::data()[4]; }

private:
    Distortion distortion_;
};

namespace util
{
template<typename T>
void unitizeZ(const T* xyz, T* uv1, T* jacob = NULL)
{
    uv1[0] = xyz[0] / xyz[2];
    uv1[1] = xyz[1] / xyz[2];

    if (jacob)
    {
        jacob[0] = T(1)/xyz[2];
        jacob[1] = T(0);
        jacob[2] = T(0);
        jacob[3] = T(1)/xyz[2];
        jacob[4] = -xyz[0]/xyz[2]/xyz[2];
        jacob[5] = -xyz[1]/xyz[2]/xyz[2];
    }
}
template<typename T>
void multplyMat3x3x3(const T* a, const T* b, T* c)
{
    c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
    c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
    c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];
    c[3] = a[0]*b[3] + a[3]*b[4] + a[6]*b[5];
    c[4] = a[1]*b[3] + a[4]*b[4] + a[7]*b[5];
    c[5] = a[2]*b[3] + a[5]*b[4] + a[8]*b[5];
    c[6] = a[0]*b[6] + a[3]*b[7] + a[6]*b[8];
    c[7] = a[1]*b[6] + a[4]*b[7] + a[7]*b[8];
    c[8] = a[2]*b[6] + a[5]*b[7] + a[8]*b[8];
}
template<typename T>
void multplyMat3x2x3(const T* a, const T* b, T* c)
{
    c[0] = a[0]*b[0] + a[3]*b[1];
    c[1] = a[1]*b[0] + a[4]*b[1];
    c[2] = a[2]*b[0] + a[5]*b[1];
    c[3] = a[0]*b[2] + a[3]*b[3];
    c[4] = a[1]*b[2] + a[4]*b[3];
    c[5] = a[2]*b[2] + a[5]*b[3];
    c[6] = a[0]*b[4] + a[3]*b[5];
    c[7] = a[1]*b[4] + a[4]*b[5];
    c[8] = a[2]*b[4] + a[5]*b[5];
}
template<typename T>
void multplyMat3x3xN(const T* a, const T* b, int N, T* c)
{
    for (int i = 0; i < N; ++i)
    {
        c[i * 3 + 0] = a[0]*b[i * 3 + 0] + a[3]*b[i * 3 + 1] + a[6]*b[i * 3 + 2];
        c[i * 3 + 1] = a[1]*b[i * 3 + 0] + a[4]*b[i * 3 + 1] + a[7]*b[i * 3 + 2];
        c[i * 3 + 2] = a[2]*b[i * 3 + 0] + a[5]*b[i * 3 + 1] + a[8]*b[i * 3 + 2];
    }
}
};


template<typename T>
class DistortionScalar
{
public:
    DistortionScalar(T* data) {}
    void project(const T* incident, T* refracted, T* jacobDistort = NULL, T* jacobRay = NULL) const;
};

template<typename PrevDistortion = DistortionScalar<double> >
class BrownDistortion : public PrevDistortion
{
    typedef BrownDistortion<PrevDistortion> Self;
public:
    typedef typename Traits<Self>::Scalar T;
    enum { kDataSize = Traits<Self>::kDataSize,
        kMyDataSize = Traits<Self>::kMyDataSize };
private:
    BrownDistortion() : data_(NULL) {}
    T* data_;
public:
    explicit BrownDistortion(T* data) : data_(data), PrevDistortion(data + kMyDataSize) {}
    void project(const T* incident, T* refracted, T* jacobDistort = NULL, T* jacobRay = NULL) const;
    const T& k1() const { return data_[0]; }
    const T& k2() const { return data_[1]; }
    const T& p1() const { return data_[2]; }
    const T& p2() const { return data_[3]; }
    const T& k3() const { return data_[4]; }
    T* data() { return data_; }
    const T* data() const { return data_; }
};

template<typename PrevDistortion = DistortionScalar<double> >
class CatadioptricDistortion : public PrevDistortion
{
    typedef CatadioptricDistortion<PrevDistortion> Self;
public:
    typedef typename Traits<Self>::Scalar T;
    enum { kDataSize = Traits<Self>::kDataSize, 
        kMyDataSize = Traits<Self>::kMyDataSize };
private:
    CatadioptricDistortion() : data_(NULL) {}
    T* data_;
public:
    explicit CatadioptricDistortion(T* data) : data_(data), PrevDistortion(data + kMyDataSize) {}
    void project(const T* incident, T* refracted, T* jacobDistort = NULL, T* jacobRay = NULL) const;
    const T& xi() const { return data_[0]; }
    T* data() { return data_; }
    const T* data() const { return data_; }
};

template<typename T>
struct Traits<DistortionScalar<T> >
{
    enum { kDataSize = 0 };
    typedef T Scalar;
};

template<typename PrevDistortion>
struct Traits<BrownDistortion<PrevDistortion> >
{
    enum { kMyDataSize = 5,
        kDataSize = kMyDataSize + Traits<PrevDistortion>::kDataSize };
    typedef typename Traits<PrevDistortion>::Scalar Scalar;
};

template<typename PrevDistortion>
struct Traits<CatadioptricDistortion<PrevDistortion> >
{
    enum { kMyDataSize = 1,
        kDataSize = kMyDataSize + Traits<PrevDistortion>::kDataSize };
    typedef typename Traits<PrevDistortion>::Scalar Scalar;
};

template<typename Distortion,
    template<typename, int> class _Storage>
struct Traits<Camera<Distortion, _Storage> >
{
    enum { kMyDataSize = 5,
        kDataSize = kMyDataSize + Traits<Distortion>::kDataSize };
    typedef typename Traits<Distortion>::Scalar Scalar;
    typedef _Storage<Scalar, kDataSize> Storage;
};

template<typename T>
void DistortionScalar<T>::project(const T* incident, T* refracted, T* jacobDistort, T* jacobRay) const
{
    memcpy(refracted, incident, 3 * sizeof(T));
    if (jacobRay != NULL)
    {
        for (int i = 0; i < 9; ++i)
            if (i % 4 == 0) jacobRay[i] = T(1);
            else jacobRay[i] = T(0);
    }
}

template<typename PrevDistortion>
void BrownDistortion<PrevDistortion>::project(const T* _incident, T* refracted, T* jacobDistort, T* jacobRay) const
{
    enum { prevJacobDistortSize = 3 * Traits<PrevDistortion>::kDataSize };
    const int prevJacobRaySize = 3 * 3;
    T prevJacobDistort[prevJacobDistortSize], prevJacobRay[prevJacobRaySize];
    T incident[3];
    PrevDistortion::project(_incident, incident, jacobDistort ? prevJacobDistort : NULL, jacobRay ? prevJacobRay : NULL);

    T uv1[3], jacobUnitize[2 * 3];
    util::unitizeZ(incident, uv1, (jacobDistort != NULL | jacobRay != NULL) ? jacobUnitize : NULL);

    const T& k1 = this->k1();
    const T& k2 = this->k2();
    const T& p1 = this->p1();
    const T& p2 = this->p2();
    const T& k3 = this->k3();

    const T& u = uv1[0];
    const T& v = uv1[1];
    const T u2 = u * u;
    const T v2 = v * v;
    const T r2 = u2 + v2;
    const T r4 = r2 * r2;
    const T r6 = r2 * r4;
    const T k = T(1) + k1 * r2 + k2 * r4 + k3 * r6;
    refracted[0] = u * k + T(2) * p1 * u * v + p2 * (r2 + T(2) * u2);
    refracted[1] = v * k + p1 * (r2 + T(2) * v2) + T(2) * p2 * u * v;
    refracted[2] = T(1);

    if (jacobDistort != NULL | jacobRay != NULL)
    {
        T jacobuv[3 * 2];
        jacobuv[0] = k2*r4 + k3*r6 + T(6)*p2*u + T(2)*p1*v + u*(T(2)*k1*u + T(4)*k2*u*r2 + T(6)*k3*u*r4) + k1*r2 + T(1);
        jacobuv[1] = T(2)*p1*u + T(2)*p2*v + v*(T(2)*k1*u + T(4)*k2*u*r2 + T(6)*k3*u*r4);
        jacobuv[2] = T(0);
        jacobuv[3] = T(2)*p1*u + T(2)*p2*v + u*(T(2)*k1*v + T(4)*k2*v*r2 + T(6)*k3*v*r4);
        jacobuv[4] = k2*r4 + k3*r6 + T(2)*p2*u + T(6)*p1*v + v*(T(2)*k1*v + T(4)*k2*v*r2 + T(6)*k3*v*r4) + k1*r2 + T(1);
        jacobuv[5] = T(0);

        T jacobIncident[9];
        util::multplyMat3x2x3(jacobuv, jacobUnitize, jacobIncident);
        if (jacobRay) util::multplyMat3x3x3(jacobIncident, prevJacobRay, jacobRay);
        if (jacobDistort)
        {
            const int jacobDistortSize = 3 * 5;

            jacobDistort[0] = u*r2;
            jacobDistort[1] = v*r2;
            jacobDistort[2] = T(0);
            jacobDistort[3] = u*r4;
            jacobDistort[4] = v*r4;
            jacobDistort[5] = T(0);
            jacobDistort[6] = T(2)*u*v;
            jacobDistort[7] = u*u + T(3)*v*v;
            jacobDistort[8] = T(0);
            jacobDistort[9] = T(3)*u*u + v*v;
            jacobDistort[10] = T(2)*u*v;
            jacobDistort[11] = T(0);
            jacobDistort[12] = u*r6;
            jacobDistort[13] = v*r6;
            jacobDistort[14] = T(0);

            T jacobPrevDistort[prevJacobDistortSize];
            util::multplyMat3x3xN(jacobIncident, prevJacobDistort, prevJacobDistortSize, jacobPrevDistort);
            memcpy(jacobDistort + jacobDistortSize, jacobPrevDistort, sizeof(T) * prevJacobDistortSize);
        }
    }
}

template<typename PrevDistortion>
void CatadioptricDistortion<PrevDistortion>::project(const T* _incident, T* refracted, T* jacobDistort, T* jacobRay) const
{
    enum { prevJacobDistortSize = 3 * Traits<PrevDistortion>::kDataSize };
    const int prevJacobRaySize = 3 * 3;
    T prevJacobDistort[prevJacobDistortSize], prevJacobRay[prevJacobRaySize];
    T incident[3];
    PrevDistortion::project(_incident, incident, jacobDistort ? prevJacobDistort : NULL, jacobRay ? prevJacobRay : NULL);

    const T& x = incident[0];
    const T& y = incident[1];
    const T& z = incident[2];
    const T r2 = x * x + y * y + z * z;
    const T r = sqrt(r2);
    const T& xi = this->xi();
    const T tmp = xi + z/r;
    refracted[0] = x/(tmp*r);
    refracted[1] = y/(tmp*r);
    refracted[2] = T(1);

    if (jacobDistort != NULL | jacobRay != NULL)
    {
        const T x2 = x * x;
        const T y2 = y * y;
        const T z2 = z * z;
        const T r3 = r2 * r;
        const T r4 = r3 * r;
        const T tmp2 = tmp * tmp;
        T jacobIncident[3 * 3];
        jacobIncident[0] = T(1)/(tmp*r) - x2/(tmp*r3) + (x2*z)/(tmp2*r4);
        jacobIncident[1] = (x*y*z)/(tmp2*r4) - (x*y)/(tmp*r3);
        jacobIncident[2] = T(0);
        jacobIncident[3] = (x*y*z)/(tmp2*r4) - (x*y)/(tmp*r3);
        jacobIncident[4] = T(1)/(tmp*r) - y2/(tmp*r3) + (y2*z)/(tmp2*r4);
        jacobIncident[5] = T(0);
        jacobIncident[6] = - (x*(T(1)/r - z2/r3))/(tmp2*r) - (x*z)/(tmp*r3);
        jacobIncident[7] = - (y*(T(1)/r - z2/r3))/(tmp2*r) - (y*z)/(tmp*r3);
        jacobIncident[8] = T(0);
        if (jacobRay) util::multplyMat3x3x3(jacobIncident, prevJacobRay, jacobRay);
        if (jacobDistort)
        {
            const int jacobDistortSize = 3 * 1;
            jacobDistort[0] = -x/(tmp2*r);
            jacobDistort[1] = -y/(tmp2*r);
            jacobDistort[2] = T(0);
            T jacobPrevDistort[prevJacobDistortSize];
            util::multplyMat3x3xN(jacobIncident, prevJacobDistort, prevJacobDistortSize, jacobPrevDistort);
            memcpy(jacobDistort + jacobDistortSize, jacobPrevDistort, sizeof(T) * prevJacobDistortSize);
        }
    }

}

typedef DistortionScalar<double> DistortionScalard;
typedef BrownDistortion<DistortionScalard> BrownDistortiond;
typedef CatadioptricDistortion<DistortionScalard> CatadioptricDistortiond;
typedef BrownDistortion<CatadioptricDistortiond> BrownCatadioptricDistortiond;

} // geometry
} // idl
#endif
