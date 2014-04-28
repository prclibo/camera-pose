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

template<bool T>
class Assert;
template<>
class Assert<true> {};

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

/*
template<typename T, 
    template<typename, template<typename, int> class> class _Distortion,
    template<typename, int> class _Storage = OwnedStorage>
class Camera
{
    Camera();
    Camera(T* data);
    Camera(const T* data);
    Camera(const Self& other) : Storage(other), orientation_(Storage::data()), location_(Storage::data()) {}
};*/

namespace util
{
template<typename T>
void unitizeZ(const T* xyz, T* xy1, T* jacob = NULL)
{
    xy1[0] = xyz[0] / xyz[2];
    xy1[1] = xyz[1] / xyz[2];

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
class BrownRadialDistortion : public Traits<BrownRadialDistortion<PrevDistortion> >::Storage
{
    typedef BrownRadialDistortion<PrevDistortion> Self;
    typedef typename Traits<Self>::Scalar T;
    typedef typename Traits<Self>::Storage Storage;
    enum { kDataSize = Traits<Self>::kDataSize };
private:
    BrownRadialDistortion() : prev_(NULL) {}
    PrevDistortion prev_;
public:
    explicit BrownRadialDistortion(T* data) : Storage(data), prev_(data + kDataSize) {}
    void project(const T* incident, T* refracted, T* jacobDistort = NULL, T* jacobRay = NULL) const;
};

template<typename PrevDistortion = DistortionScalar<double> >
class BrownDistortion : public Traits<BrownDistortion<PrevDistortion> >::Storage
{
    typedef BrownDistortion<PrevDistortion> Self;
    typedef typename Traits<Self>::Scalar T;
    typedef typename Traits<Self>::Storage Storage;
    enum { kDataSize = Traits<Self>::kDataSize };
private:
    BrownDistortion() : prev_(NULL), radial_(NULL) {}
    BrownRadialDistortion<PrevDistortion> radial_;
    PrevDistortion prev_;
public:
    explicit BrownDistortion(T* data) : Storage(data), prev_(data + kDataSize), radial_(data) {}
    void project(const T* incident, T* refracted, T* jacobDistort = NULL, T* jacobRay = NULL) const;
};

template<typename T>
struct Traits<DistortionScalar<T> >
{
    enum { kDataSize = 0 };
    typedef T Scalar;
};

template<typename PrevDistortion>
struct Traits<BrownRadialDistortion<PrevDistortion> >
{
    enum { kMyDataSize = 3,
        kDataSize = kMyDataSize + Traits<PrevDistortion>::kDataSize };
    typedef typename Traits<PrevDistortion>::Scalar Scalar;
    typedef SharedStorage<Scalar, kDataSize> Storage;
};

template<typename PrevDistortion>
struct Traits<BrownDistortion<PrevDistortion> >
{
    enum { kMyDataSize = 2,
        kDataSize = kMyDataSize + Traits<BrownRadialDistortion<PrevDistortion> >::kDataSize };
    typedef typename Traits<PrevDistortion>::Scalar Scalar;
    typedef SharedStorage<Scalar, kDataSize> Storage;
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
void BrownRadialDistortion<PrevDistortion>::project(const T* _incident, T* refracted, T* jacobDistort, T* jacobRay) const
{
    enum { prevJacobDistortSize = 3 * Traits<PrevDistortion>::kDataSize };
    const int prevJacobRaySize = 3 * 3;
    T prevJacobDistort[prevJacobDistortSize], prevJacobRay[prevJacobRaySize];
    T incident[3];
    prev_.project(_incident, incident, jacobDistort ? prevJacobDistort : NULL, jacobRay ? prevJacobRay : NULL);

    T xy1[3], jacobUnitize[2 * 3];
    util::unitizeZ(incident, xy1, (jacobDistort != NULL | jacobRay != NULL) ? jacobUnitize : NULL);

    const T& k1 = Storage::data()[0];
    const T& k2 = Storage::data()[1];
    const T& k3 = Storage::data()[2];

    const T& x = xy1[0];
    const T& y = xy1[1];
    const T r2 = x * x + y * y;
    const T r4 = r2 * r2;
    const T r6 = r2 * r4;
    const T k = T(1) + k1 * r2 + k2 * r4 + k3 * r6;
    refracted[0] = x * k;
    refracted[1] = y * k;
    refracted[2] = T(1);

    if (jacobDistort != NULL | jacobRay != NULL)
    {
        T jacobXy[3 * 2];
        jacobXy[0] = k2*r4 + k3*r6 + x*(T(2)*k1*x + T(4)*k2*x*r2 + T(6)*k3*x*r4) + k1*r2 + T(1);
        jacobXy[1] = y*(T(2)*k1*x + T(4)*k2*x*r2 + T(6)*k3*x*r4);
        jacobXy[2] = T(0);
        jacobXy[3] = x*(T(2)*k1*y + T(4)*k2*y*r2 + T(6)*k3*y*r4);
        jacobXy[4] = k2*r4 + k3*r6 + y*(T(2)*k1*y + T(4)*k2*y*r2 + T(6)*k3*y*r4) + k1*r2 + T(1);
        jacobXy[5] = T(0);

        T jacobIncident[9];
        util::multplyMat3x2x3(jacobXy, jacobUnitize, jacobIncident);
        if (jacobRay) util::multplyMat3x3x3(jacobIncident, prevJacobRay, jacobRay);
        if (jacobDistort)
        {
            T jacobXyDistort[6];
            jacobXyDistort[0] = x*r2;
            jacobXyDistort[1] = y*r4;
            jacobXyDistort[2] = x*r6;
            jacobXyDistort[3] = y*r2;
            jacobXyDistort[4] = x*r4;
            jacobXyDistort[5] = y*r6;
            util::multplyMat3x2x3(jacobXy, jacobXyDistort, jacobDistort);

            T jacobPrevDistort[prevJacobDistortSize];
            util::multplyMat3x3xN(jacobIncident, prevJacobDistort, prevJacobDistortSize, jacobPrevDistort);
            memcpy(jacobDistort + 9, jacobDistort, sizeof(T) * prevJacobDistortSize);
        }
    }
}

template<typename PrevDistortion>
void BrownDistortion<PrevDistortion>::project(const T* _incident, T* refracted, T* jacobDistort, T* jacobRay) const
{
    const int prevJacobDistortSize = 3 * Traits<PrevDistortion>::kDataSize;
    T incident[3];
    prev_.project(_incident, incident, jacobDistort ? jacobDistort : NULL, jacobRay ? jacobRay : NULL);

    const int radialJacobDistortSize = 3 * Traits<BrownRadialDistortion<PrevDistortion> >::kMyDataSize;
    T radialRefracted[3];

    radial_.project(incident, radialRefracted, jacobDistort ? jacobDistort : NULL, 
            jacobRay ? jacobRay : NULL);
    assert(radialRefracted[2] == T(1));

    T xy1[3], jacobUnitize[2 * 3];
    util::unitizeZ(incident, xy1, (jacobDistort != NULL | jacobRay != NULL) ? jacobUnitize : NULL);

    const T& p1 = Storage::data()[0];
    const T& p2 = Storage::data()[1];
    T x = xy1[0];
    T y = xy1[1];
    T r2 = x * x + y * y;

    T dx = T(2) * p1 * x * y + p2 * (r2 + T(2) * x * x);
    T dy = p1 * (r2 + T(2) * y * y) + T(2) * p2 * x * y;

    refracted[0] = radialRefracted[0] / radialRefracted[2] + dx;
    refracted[1] = radialRefracted[1] / radialRefracted[2] + dy;
    refracted[2] = T(1);

    T* jacobTangent = jacobDistort + radialJacobDistortSize;


}

typedef DistortionScalar<double> DistortionScalard;
typedef BrownRadialDistortion<DistortionScalard> BrownRadialDistortiond;
typedef BrownDistortion<DistortionScalard> BrownDistortiond;
typedef BrownDistortion<CatadioptricDistortion<DistortionScalar<double> > >

} // geometry
} // idl
#endif
