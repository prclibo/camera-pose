// Ask libo for details.

#ifndef POSE_H_
#define POSE_H_

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

template<typename Derived, template<typename, int> class Storage> class AngleAxis;
template<typename Derived, template<typename, int> class Storage> class Quaternion;
template<typename Derived, template<typename, int> class Storage> class RotationMatrix;

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

template<typename Derived>
class OrientationBase : public Traits<Derived>::Storage
{
    typedef typename Traits<Derived>::Scalar T;
    typedef typename Traits<Derived>::Storage Storage;
    typedef typename Traits<Derived>::OwnedType OwnedDerived;
public:
    enum {
        kDataSize = Traits<Derived>::kDataSize,
    };

    OrientationBase() : Storage() {}
    OrientationBase(T* data) : Storage(data) {}
    OrientationBase(const T* data) : Storage(data) {}
    
    // Self& operator =(const Self& other) { Storage::operator =(other); return *this; }

    inline const Derived& derived() const { return *static_cast<const Derived*>(this); }
    inline Derived& derived() { return *static_cast<Derived*>(this); }
    inline static OwnedDerived Identity() { return Derived::Identity(); }
    inline void transformIn(const T* point, T* transformed, T* jacob = NULL) const;
    inline void transformOut(const T* point, T* transformed, T* jacob = NULL) const;
    inline OwnedDerived combine(const Derived& other) const;
    inline OwnedDerived relate(const Derived& other) const;
    inline OwnedDerived inverse() const { derived().inverse(); }
    inline void toRotationMatrix(T* matrix, T* jacob = NULL) const { derived().toRotationMatrix(); }
protected:
    ~OrientationBase() {}
};

template<typename T, template<typename, int> class _Storage = OwnedStorage>
class RotationMatrix : public OrientationBase<RotationMatrix<T, _Storage> >
{
public:
    typedef RotationMatrix<T, _Storage> Self;
    typedef OrientationBase<Self> BaseType;
    typedef typename Traits<Self>::Storage Storage;
    typedef typename Traits<Self>::OwnedType OwnedSelf;
public:
    RotationMatrix();
    RotationMatrix(T* data) : BaseType(data) {}
    RotationMatrix(const T* data) : BaseType(data) {}
    // RotationMatrix(const Self& other) { CheckOwnedStorage<_Storage> this_constructor_is_only_for_OwnedStoraged_based_object; *this = other; }
    template<typename Other>
    RotationMatrix(const Other& other) { CheckOwnedStorage<_Storage> this_constructor_is_only_for_OwnedStoraged_based_object; *this = other; }

    template<typename NewT>
    RotationMatrix<NewT> cast() const;
    // Self& operator =(const Self& other) { BaseType::operator =(other); return *this; };
    template<typename Other>
    Self& operator =(const Other& other);
    
    static OwnedSelf Identity();
    OwnedSelf inverse() const;
    void toRotationMatrix(T* matrix/*[3 * 3]*/, T* jacob/*[9 * 9]*/ = NULL) const;
    void toAngleAxis(T* angleAxis/*[3]*/, T* jacob/*[3 * 9]*/ = NULL) const;
};

template<typename T, template<typename, int> class _Storage = OwnedStorage>
class AngleAxis : public OrientationBase<AngleAxis<T, _Storage> >
{
public:
    typedef AngleAxis<T, _Storage> Self;
    typedef OrientationBase<Self> BaseType;
    typedef typename Traits<Self>::Storage Storage;
    typedef typename Traits<Self>::OwnedType OwnedSelf;
public:
    AngleAxis();
    AngleAxis(T* data) : BaseType(data) {}
    AngleAxis(const T* data) : BaseType(data) {}
    // AngleAxis(const Self& other) { CheckOwnedStorage<_Storage> this_constructor_is_only_for_OwnedStoraged_based_object; *this = other; }
    template<typename Other>
    AngleAxis(const Other& other) { CheckOwnedStorage<_Storage> this_constructor_is_only_for_OwnedStoraged_based_object; *this = other; }
    AngleAxis(T x, T y, T z);

    template<typename NewT>
    AngleAxis<NewT> cast() const;
    // Self& operator =(const Self& other) { BaseType::operator =(other); return *this; };
    template<typename Other>
    Self& operator =(const Other& other);

    static OwnedSelf Identity();
    OwnedSelf inverse() const;
    void toRotationMatrix(T* rotationMatrix/*[3 * 3]*/, T* jacob/*[3 * 9]*/ = NULL) const;
    void toAngleAxis(T* angleAxis/*[3]*/, T* jacob/*[3 * 9]*/ = NULL) const;
};

template<typename Derived>
class PoseBase : public Traits<Derived>::Storage
{
    typedef typename Traits<Derived>::Orientation Orientation;
    typedef typename Traits<Derived>::Storage Storage;
    typedef typename Traits<Derived>::Scalar T;
    typedef typename Traits<Derived>::OwnedType OwnedDerived;
    enum { 
        kOrientSize = Traits<Orientation>::kDataSize,
        kDataSize = kOrientSize + 3
    };
    typedef PoseBase<Derived> Self;
protected:
    Orientation orientation_;
    T* location_;

public:
    PoseBase();
    PoseBase(T* data);
    PoseBase(const T* data);
    PoseBase(const Self& other) : Storage(other), orientation_(Storage::data()), location_(Storage::data()) {}
    
    Orientation& orientation() { return orientation_; }
    const Orientation& orientation() const { return orientation_; }
    T* location() { return Storage::data() + kOrientSize; }
    const T* location() const { return Storage::data() + kOrientSize; }

    inline const Derived& derived() const { return *static_cast<const Derived*>(this); }
    inline Derived& derived() { return *static_cast<Derived*>(this); }
    inline void transformIn(const T* point, T* transformed, T* jacob = NULL) const;
    inline void transformOut(const T* point, T* transformed, T* jacob = NULL) const;
    inline OwnedDerived combine(const Derived& other) const;
    inline OwnedDerived relate(const Derived& other) const;
    inline OwnedDerived inverse() const;
    void toMatrix3x4(T* matrix/*[4 * 4]*/, T* jacob/*[16 * kDataSize]*/ = NULL) const;
    
    static OwnedDerived Identity();
};

template<typename T,
    template<typename, template<typename, int> class> class _Orientation,
    template<typename, int> class _Storage = OwnedStorage>
class Pose : public PoseBase<Pose<T, _Orientation, _Storage> > 
{
public:
    typedef _Orientation<T, _Storage> Orientation;
    enum { 
        kOrientSize = Traits<Orientation>::kDataSize,
        kDataSize = kOrientSize + 3
    };
    typedef _Storage<T, kDataSize> Storage;
    typedef Pose<T, _Orientation, _Storage> Self;
    typedef PoseBase<Self> BaseType;
public:
    Pose() : BaseType() {}
    Pose(T* data) : BaseType(data) {}
    Pose(const T* data) : BaseType(data) {}
    // Pose(const Self& other) : BaseType(other) { CheckOwnedStorage<_Storage> this_constructor_is_only_for_OwnedStoraged_based_object; }
    Pose(const Orientation& orient, const T* loc);

    template<typename Other>
    Pose(const Other& other) { CheckOwnedStorage<_Storage> this_constructor_is_only_for_OwnedStoraged_based_object; *this = other; }

    template<typename NewT>
    Pose<NewT, _Orientation> cast() const;
    template<typename Other>
    Self& operator =(const Other& other);
};

template<typename T, template<typename, int> class _Storage>
struct Traits<RotationMatrix<T, _Storage > >
{
    enum
    {
        kDataSize = 9,
    };
    typedef T Scalar;
    typedef _Storage<T, 9> Storage;
    typedef RotationMatrix<T, OwnedStorage> OwnedType;
};

template<typename T, template<typename, int> class _Storage>
struct Traits<AngleAxis<T, _Storage > >
{
    enum
    {
        kDataSize = 3,
    };
    typedef T Scalar;
    typedef _Storage<T, 3> Storage;
    typedef AngleAxis<T, OwnedStorage> OwnedType;
};

template<typename T,
    template<typename, template<typename, int> class> class _Orientation,
    template<typename, int> class _Storage>
struct Traits<Pose<T, _Orientation, _Storage> >
{
    typedef _Orientation<T, SharedStorage> Orientation;
    enum { 
        kOrientSize = Traits<Orientation>::kDataSize,
        kDataSize = kOrientSize + 3
    };
    typedef T Scalar;
    typedef _Storage<T, kDataSize> Storage;
    typedef Pose<T, _Orientation, OwnedStorage> OwnedType;
};
// ---------------------------------------------------------------------------

template<typename Derived>
void OrientationBase<Derived>::transformOut(const T* point, T* transformed, T* jacob) const
{
    T H[3 * 3];
    derived().toRotationMatrix(H);
    transformed[0] = H[0]*point[0] + H[3]*point[1] + H[6]*point[2];
    transformed[1] = H[1]*point[0] + H[4]*point[1] + H[7]*point[2];
    transformed[2] = H[2]*point[0] + H[5]*point[1] + H[8]*point[2];
    if (jacob)
    {
        jacob[0] = point[0];
        jacob[1] = T(0);
        jacob[2] = T(0);
        jacob[3] = T(0);
        jacob[4] = point[0];
        jacob[5] = T(0);
        jacob[6] = T(0);
        jacob[7] = T(0);
        jacob[8] = point[0];
        jacob[9] = point[1];
        jacob[10] = T(0);
        jacob[11] = T(0);
        jacob[12] = T(0);
        jacob[13] = point[1];
        jacob[14] = T(0);
        jacob[15] = T(0);
        jacob[16] = T(0);
        jacob[17] = point[1];
        jacob[18] = point[2];
        jacob[19] = T(0);
        jacob[20] = T(0);
        jacob[21] = T(0);
        jacob[22] = point[2];
        jacob[23] = T(0);
        jacob[24] = T(0);
        jacob[25] = T(0);
        jacob[26] = point[2];
    }
}

template<typename Derived>
void OrientationBase<Derived>::transformIn(const T* point, T* transformed, T* jacob) const
{
    derived().inverse().transformOut(point, transformed, jacob);
}

template<typename Derived>
typename OrientationBase<Derived>::OwnedDerived OrientationBase<Derived>::combine(const Derived& other) const
{
    T H0[3 * 3], H1[3 * 3], C[3 * 3];
    derived().toRotationMatrix(H0);
    other.toRotationMatrix(H1);

    C[0] = H0[0]*H1[0] + H1[1]*H0[3] + H1[2]*H0[6];
    C[1] = H0[1]*H1[0] + H1[1]*H0[4] + H1[2]*H0[7];
    C[2] = H1[0]*H0[2] + H1[1]*H0[5] + H1[2]*H0[8];
    C[3] = H0[0]*H1[3] + H0[3]*H1[4] + H0[6]*H1[5];
    C[4] = H0[1]*H1[3] + H0[4]*H1[4] + H1[5]*H0[7];
    C[5] = H0[2]*H1[3] + H0[5]*H1[4] + H1[5]*H0[8];
    C[6] = H0[0]*H1[6] + H0[3]*H1[7] + H0[6]*H1[8];
    C[7] = H0[1]*H1[6] + H0[4]*H1[7] + H0[7]*H1[8];
    C[8] = H0[2]*H1[6] + H0[5]*H1[7] + H0[8]*H1[8];

    RotationMatrix<T>  matrix(C);
    return Derived(matrix);
}

template<typename Derived>
typename OrientationBase<Derived>::OwnedDerived OrientationBase<Derived>::relate(const Derived& other) const
{
    T H0[3 * 3], H1[3 * 3], C[3 * 3];
    derived().inverse().toRotationMatrix(H0);
    other.toRotationMatrix(H1);

    C[0] = H0[0]*H1[0] + H1[1]*H0[3] + H1[2]*H0[6];
    C[1] = H0[1]*H1[0] + H1[1]*H0[4] + H1[2]*H0[7];
    C[2] = H1[0]*H0[2] + H1[1]*H0[5] + H1[2]*H0[8];
    C[3] = H0[0]*H1[3] + H0[3]*H1[4] + H0[6]*H1[5];
    C[4] = H0[1]*H1[3] + H0[4]*H1[4] + H1[5]*H0[7];
    C[5] = H0[2]*H1[3] + H0[5]*H1[4] + H1[5]*H0[8];
    C[6] = H0[0]*H1[6] + H0[3]*H1[7] + H0[6]*H1[8];
    C[7] = H0[1]*H1[6] + H0[4]*H1[7] + H0[7]*H1[8];
    C[8] = H0[2]*H1[6] + H0[5]*H1[7] + H0[8]*H1[8];

    RotationMatrix<T> matrix(C);
    return Derived(matrix);
}

template<typename T, template<typename, int> class _Storage>
RotationMatrix<T, _Storage>::RotationMatrix()
: OrientationBase<RotationMatrix<T, _Storage> >::OrientationBase()
{
    OrientationBase<RotationMatrix<T, _Storage> >::data()[0] = T(1);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[1] = T(0);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[2] = T(0);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[3] = T(0);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[4] = T(1);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[5] = T(0);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[6] = T(0);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[7] = T(0);
    OrientationBase<RotationMatrix<T, _Storage> >::data()[8] = T(1);
}

template<typename T, template<typename, int> class _Storage>
template<typename NewT>
RotationMatrix<NewT> RotationMatrix<T, _Storage>::cast() const
{
    RotationMatrix<NewT> ret;
    for (int i = 0; i < Traits<Self>::kDataSize; ++i)
        ret.data()[i] = NewT(Storage::data()[i]);
    return ret;
}

template<typename T, template<typename, int> class _Storage>
template<typename Other>
RotationMatrix<T, _Storage>& RotationMatrix<T, _Storage>::operator =(const Other& other)
{
    T m[9];
    other.template cast<T>().toRotationMatrix(m);
    memcpy(Storage::data(), m, Traits<Self>::kDataSize * sizeof(T));
    return *this;
}

template<typename T, template<typename, int> class _Storage>
void RotationMatrix<T, _Storage>::toRotationMatrix(T* matrix, T* jacob) const
{
    memcpy(matrix, OrientationBase<RotationMatrix<T, _Storage> >::data(),
            sizeof(T) * OrientationBase<RotationMatrix<T, _Storage> >::dataSize());

    if (!jacob) return;
    for (int i = 0; i < 9 * 9; ++i)
        jacob[i] = (i % 10 == 0) ? T(1) : T(0);
}

template<typename T, template<typename, int> class _Storage>
void RotationMatrix<T, _Storage>::toAngleAxis(T* angleAxis, T* jacob) const
{
    const T* R = Storage::data();
    T* aa = angleAxis;

    T r13 = (R[1]/T(2) - R[3]/T(2));
    T r26 = (R[2]/T(2) - R[6]/T(2));
    T r57 = (R[5]/T(2) - R[7]/T(2));
    T sint2 = r13 * r13 + r26 * r26 + r57 * r57;
    if (sint2 > T(1))  sint2 = T(1);
    T sint = sqrt(sint2);
    T sint3 = sint2 * sint;
    T cost = (R[0] + R[4] + R[8] - T(1)) / T(2);
    if (cost > T(1)) cost = T(1);
    if (cost < T(-1)) cost = T(-1);
    T abscost = cost > T(0) ? cost : -cost;
    T abssint = sint > T(0) ? sint : -sint;
    T cost2 = cost * cost;

    T theta = cost > 0 ? asin(sint) : T(M_PI) - asin(sint);

    aa[0] = r57 / sint * theta;
    aa[1] = -r26 / sint * theta;
    aa[2] = r13 / sint * theta;

    if (jacob && cost >= 0)
    {
        jacob[0] = 0;
        jacob[1] = 0;
        jacob[2] = 0;
        jacob[3] = (r13*r57)/(2*sint2*abscost) - (theta*r13*r57)/(2*sint3);
        jacob[4] = (theta*r13*r26)/(2*sint3) - (r13*r26)/(2*sint2*abscost);
        jacob[5] = theta/(2*sint) + r13*r13/(2*sint2*abscost) - (theta*r13*r13)/(2*sint3);
        jacob[6] = (r26*r57)/(2*sint2*abscost) - (theta*r26*r57)/(2*sint3);
        jacob[7] = (theta*r26*r26)/(2*sint3) - r26*r26/(2*sint2*abscost) - theta/(2*sint);
        jacob[8] = (r13*r26)/(2*sint2*abscost) - (theta*r13*r26)/(2*sint3);
        jacob[9] = (theta*r13*r57)/(2*sint3) - (r13*r57)/(2*sint2*abscost);
        jacob[10] = (r13*r26)/(2*sint2*abscost) - (theta*r13*r26)/(2*sint3);
        jacob[11] = (theta*r13*r13)/(2*sint3) - r13*r13/(2*sint2*abscost) - theta/(2*sint);
        jacob[12] = 0;
        jacob[13] = 0;
        jacob[14] = 0;
        jacob[15] = theta/(2*sint) + r57*r57/(2*sint2*abscost) - (theta*r57*r57)/(2*sint3);
        jacob[16] = (theta*r26*r57)/(2*sint3) - (r26*r57)/(2*sint2*abscost);
        jacob[17] = (r13*r57)/(2*sint2*abscost) - (theta*r13*r57)/(2*sint3);
        jacob[18] = (theta*r26*r57)/(2*sint3) - (r26*r57)/(2*sint2*abscost);
        jacob[19] = theta/(2*sint) + r26*r26/(2*sint2*abscost) - (theta*r26*r26)/(2*sint3);
        jacob[20] = (theta*r13*r26)/(2*sint3) - (r13*r26)/(2*sint2*abscost);
        jacob[21] = (theta*r57*r57)/(2*sint3) - r57*r57/(2*sint2*abscost) - theta/(2*sint);
        jacob[22] = (r26*r57)/(2*sint2*abscost) - (theta*r26*r57)/(2*sint3);
        jacob[23] = (theta*r13*r57)/(2*sint3) - (r13*r57)/(2*sint2*abscost);
        jacob[24] = 0;
        jacob[25] = 0;
        jacob[26] = 0;
    }
    else if (jacob && cost < 0)
    {
        jacob[0] = 0;
        jacob[1] = 0;
        jacob[2] = 0;
        jacob[3] = - (r13*r57)/(2*sint2*abscost) - (r13*r57*theta)/(2*sint3);
        jacob[4] = (r13*r26)/(2*sint2*abscost) + (r13*r26*theta)/(2*sint3);
        jacob[5] = theta/(2*sint) - r13*r13/(2*sint2*abscost) - (r13*r13*theta)/(2*sint3);
        jacob[6] = - (r26*r57)/(2*sint2*abscost) - (r26*r57*theta)/(2*sint3);
        jacob[7] = r26*r26/(2*sint2*abscost) - theta/(2*sint) + (r26*r26*theta)/(2*sint3);
        jacob[8] = - (r13*r26)/(2*sint2*abscost) - (r13*r26*theta)/(2*sint3);
        jacob[9] = (r13*r57)/(2*sint2*abscost) + (r13*r57*theta)/(2*sint3);
        jacob[10] = - (r13*r26)/(2*sint2*abscost) - (r13*r26*theta)/(2*sint3);
        jacob[11] = r13*r13/(2*sint2*abscost) - theta/(2*sint) + (r13*r13*theta)/(2*sint3);
        jacob[12] = 0;
        jacob[13] = 0;
        jacob[14] = 0;
        jacob[15] = theta/(2*sint) - r57*r57/(2*sint2*abscost) - (r57*r57*theta)/(2*sint3);
        jacob[16] = (r26*r57)/(2*sint2*abscost) + (r26*r57*theta)/(2*sint3);
        jacob[17] = - (r13*r57)/(2*sint2*abscost) - (r13*r57*theta)/(2*sint3);
        jacob[18] = (r26*r57)/(2*sint2*abscost) + (r26*r57*theta)/(2*sint3);
        jacob[19] = theta/(2*sint) - r26*r26/(2*sint2*abscost) - (r26*r26*theta)/(2*sint3);
        jacob[20] = (r13*r26)/(2*sint2*abscost) + (r13*r26*theta)/(2*sint3);
        jacob[21] = r57*r57/(2*sint2*abscost) - theta/(2*sint) + (r57*r57*theta)/(2*sint3);
        jacob[22] = - (r26*r57)/(2*sint2*abscost) - (r26*r57*theta)/(2*sint3);
        jacob[23] = (r13*r57)/(2*sint2*abscost) + (r13*r57*theta)/(2*sint3);
        jacob[24] = 0;
        jacob[25] = 0;
        jacob[26] = 0;
    
    }
}

template<typename T, template<typename, int> class _Storage>
typename RotationMatrix<T, _Storage>::OwnedSelf RotationMatrix<T, _Storage>::inverse() const
{
    RotationMatrix<T, OwnedStorage> inversed(*this);
    std::swap(inversed.data()[1], inversed.data()[3]);
    std::swap(inversed.data()[2], inversed.data()[6]);
    std::swap(inversed.data()[5], inversed.data()[7]);
    return inversed;
}

template<typename T, template<typename, int> class _Storage>
typename RotationMatrix<T, _Storage>::OwnedSelf RotationMatrix<T, _Storage>::Identity()
{
    OwnedSelf ret;
    for (int i = 0; i < 9; ++i)
        ret.data()[i] = (i % 4 == 0) ? T(1):T(0);
    return ret;
}
template<typename T, template<typename, int> class _Storage>
AngleAxis<T, _Storage>::AngleAxis()
: BaseType::OrientationBase()
{
    Storage::data()[0] = T(0);
    Storage::data()[1] = T(0);
    Storage::data()[2] = T(0);
}

template<typename T, template<typename, int> class _Storage>
AngleAxis<T, _Storage>::AngleAxis(T x, T y, T z)
: BaseType::OrientationBase()
{
    Storage::data()[0] = x;
    Storage::data()[1] = y;
    Storage::data()[2] = z;
}

template<typename T, template<typename, int> class _Storage>
template<typename NewT>
AngleAxis<NewT> AngleAxis<T, _Storage>::cast() const
{
    AngleAxis<NewT> ret;
    for (int i = 0; i < Traits<Self>::kDataSize; ++i)
        ret.data()[i] = NewT(Storage::data()[i]);
    return ret;
}

template<typename T, template<typename, int> class _Storage>
void AngleAxis<T, _Storage>::toRotationMatrix(T* matrix, T* jacob) const
{
    T* R = matrix;
    const T* aa = Storage::data();
    T aa0aa0 = aa[0] * aa[0];
    T aa1aa1 = aa[1] * aa[1];
    T aa2aa2 = aa[2] * aa[2];
    T aa0aa1 = aa[0] * aa[1];
    T aa0aa2 = aa[0] * aa[2];
    T aa1aa2 = aa[1] * aa[2];
    T theta2 = aa0aa0 + aa1aa1 + aa2aa2;
    T theta = sqrt(theta2);
    T sint = sin(theta);
    T cost = cos(theta);
    static const T kOne = T(1.0);
    if (theta2 > T(std::numeric_limits<double>::epsilon())) {
        R[0] = aa0aa0/theta2 - cost*(aa0aa0/theta2 - T(1));
        R[1] = (aa0aa1)/theta2 + (sint*aa[2])/theta - (cost*aa0aa1)/theta2;
        R[2] = (aa0aa2)/theta2 - (sint*aa[1])/theta - (cost*aa0aa2)/theta2;
        R[3] = (aa0aa1)/theta2 - (sint*aa[2])/theta - (cost*aa0aa1)/theta2;
        R[4] = aa1aa1/theta2 - cost*(aa1aa1/theta2 - T(1));
        R[5] = (aa1aa2)/theta2 + (sint*aa[0])/theta - (cost*aa1aa2)/theta2;
        R[6] = (aa0aa2)/theta2 + (sint*aa[1])/theta - (cost*aa0aa2)/theta2;
        R[7] = (aa1aa2)/theta2 - (sint*aa[0])/theta - (cost*aa1aa2)/theta2;
        R[8] = aa2aa2/theta2 - cost*(aa2aa2/theta2 - T(1));
    }
    if (jacob)
    {
        // TODO(libo): More substitution of aa[?]*aa[?]*a[?]
        T theta3 = theta2 * theta;
        T theta4 = theta3 * theta;
        T aa0aa0aa0 = aa0aa0 * aa[0];
        T aa1aa1aa1 = aa1aa1 * aa[1];
        T aa2aa2aa2 = aa2aa2 * aa[2];
        jacob[0] = (2*aa[0])/theta2 - cost*((2*aa[0])/theta2 - (2*aa0aa0aa0)/theta4) - (2*aa0aa0aa0)/theta4 + (sint*aa[0]*(aa0aa0/theta2 - T(1)))/theta;
        jacob[1] = aa[1]/theta2 - (2*aa0aa0*aa[1])/theta4 - (cost*aa[1])/theta2 + (cost*aa0aa2)/theta2 - (sint*aa0aa2)/theta3 + (2*cost*aa0aa0*aa[1])/theta4 + (sint*aa0aa0*aa[1])/theta3;
        jacob[2] = aa[2]/theta2 - (2*aa0aa0*aa[2])/theta4 - (cost*aa[2])/theta2 - (cost*aa0aa1)/theta2 + (sint*aa0aa1)/theta3 + (2*cost*aa0aa0*aa[2])/theta4 + (sint*aa0aa0*aa[2])/theta3;
        jacob[3] = aa[1]/theta2 - (2*aa0aa0*aa[1])/theta4 - (cost*aa[1])/theta2 - (cost*aa0aa2)/theta2 + (sint*aa0aa2)/theta3 + (2*cost*aa0aa0*aa[1])/theta4 + (sint*aa0aa0*aa[1])/theta3;
        jacob[4] = (sint*aa[0]*(aa1aa1/theta2 - T(1)))/theta - (2*aa[0]*aa1aa1)/theta4 + (2*cost*aa[0]*aa1aa1)/theta4;
        jacob[5] = sint/theta + (cost*aa0aa0)/theta2 - (sint*aa0aa0)/theta3 - (2*aa0aa1*aa[2])/theta4 + (2*cost*aa0aa1*aa[2])/theta4 + (sint*aa0aa1*aa[2])/theta3;
        jacob[6] = aa[2]/theta2 - (2*aa0aa0*aa[2])/theta4 - (cost*aa[2])/theta2 + (cost*aa0aa1)/theta2 - (sint*aa0aa1)/theta3 + (2*cost*aa0aa0*aa[2])/theta4 + (sint*aa0aa0*aa[2])/theta3;
        jacob[7] = (sint*aa0aa0)/theta3 - (cost*aa0aa0)/theta2 - sint/theta - (2*aa0aa1*aa[2])/theta4 + (2*cost*aa0aa1*aa[2])/theta4 + (sint*aa0aa1*aa[2])/theta3;
        jacob[8] = (sint*aa[0]*(aa2aa2/theta2 - T(1)))/theta - (2*aa[0]*aa2aa2)/theta4 + (2*cost*aa[0]*aa2aa2)/theta4;
        jacob[9] = (sint*aa[1]*(aa0aa0/theta2 - T(1)))/theta - (2*aa0aa0*aa[1])/theta4 + (2*cost*aa0aa0*aa[1])/theta4;
        jacob[10] = aa[0]/theta2 - (2*aa[0]*aa1aa1)/theta4 - (cost*aa[0])/theta2 + (cost*aa1aa2)/theta2 - (sint*aa1aa2)/theta3 + (2*cost*aa[0]*aa1aa1)/theta4 + (sint*aa[0]*aa1aa1)/theta3;
        jacob[11] = (sint*aa1aa1)/theta3 - (cost*aa1aa1)/theta2 - sint/theta - (2*aa0aa1*aa[2])/theta4 + (2*cost*aa0aa1*aa[2])/theta4 + (sint*aa0aa1*aa[2])/theta3;
        jacob[12] = aa[0]/theta2 - (2*aa[0]*aa1aa1)/theta4 - (cost*aa[0])/theta2 - (cost*aa1aa2)/theta2 + (sint*aa1aa2)/theta3 + (2*cost*aa[0]*aa1aa1)/theta4 + (sint*aa[0]*aa1aa1)/theta3;
        jacob[13] = (2*aa[1])/theta2 - cost*((2*aa[1])/theta2 - (2*aa1aa1aa1)/theta4) - (2*aa1aa1aa1)/theta4 + (sint*aa[1]*(aa1aa1/theta2 - T(1)))/theta;
        jacob[14] = aa[2]/theta2 - (2*aa1aa1*aa[2])/theta4 - (cost*aa[2])/theta2 + (cost*aa0aa1)/theta2 - (sint*aa0aa1)/theta3 + (2*cost*aa1aa1*aa[2])/theta4 + (sint*aa1aa1*aa[2])/theta3;
        jacob[15] = sint/theta + (cost*aa1aa1)/theta2 - (sint*aa1aa1)/theta3 - (2*aa0aa1*aa[2])/theta4 + (2*cost*aa0aa1*aa[2])/theta4 + (sint*aa0aa1*aa[2])/theta3;
        jacob[16] = aa[2]/theta2 - (2*aa1aa1*aa[2])/theta4 - (cost*aa[2])/theta2 - (cost*aa0aa1)/theta2 + (sint*aa0aa1)/theta3 + (2*cost*aa1aa1*aa[2])/theta4 + (sint*aa1aa1*aa[2])/theta3;
        jacob[17] = (sint*aa[1]*(aa2aa2/theta2 - T(1)))/theta - (2*aa[1]*aa2aa2)/theta4 + (2*cost*aa[1]*aa2aa2)/theta4;
        jacob[18] = (sint*aa[2]*(aa0aa0/theta2 - T(1)))/theta - (2*aa0aa0*aa[2])/theta4 + (2*cost*aa0aa0*aa[2])/theta4;
        jacob[19] = sint/theta + (cost*aa2aa2)/theta2 - (sint*aa2aa2)/theta3 - (2*aa0aa1*aa[2])/theta4 + (2*cost*aa0aa1*aa[2])/theta4 + (sint*aa0aa1*aa[2])/theta3;
        jacob[20] = aa[0]/theta2 - (2*aa[0]*aa2aa2)/theta4 - (cost*aa[0])/theta2 - (cost*aa1aa2)/theta2 + (sint*aa1aa2)/theta3 + (2*cost*aa[0]*aa2aa2)/theta4 + (sint*aa[0]*aa2aa2)/theta3;
        jacob[21] = (sint*aa2aa2)/theta3 - (cost*aa2aa2)/theta2 - sint/theta - (2*aa0aa1*aa[2])/theta4 + (2*cost*aa0aa1*aa[2])/theta4 + (sint*aa0aa1*aa[2])/theta3;
        jacob[22] = (sint*aa[2]*(aa1aa1/theta2 - T(1)))/theta - (2*aa1aa1*aa[2])/theta4 + (2*cost*aa1aa1*aa[2])/theta4;
        jacob[23] = aa[1]/theta2 - (2*aa[1]*aa2aa2)/theta4 - (cost*aa[1])/theta2 + (cost*aa0aa2)/theta2 - (sint*aa0aa2)/theta3 + (2*cost*aa[1]*aa2aa2)/theta4 + (sint*aa[1]*aa2aa2)/theta3;
        jacob[24] = aa[0]/theta2 - (2*aa[0]*aa2aa2)/theta4 - (cost*aa[0])/theta2 + (cost*aa1aa2)/theta2 - (sint*aa1aa2)/theta3 + (2*cost*aa[0]*aa2aa2)/theta4 + (sint*aa[0]*aa2aa2)/theta3;
        jacob[25] = aa[1]/theta2 - (2*aa[1]*aa2aa2)/theta4 - (cost*aa[1])/theta2 - (cost*aa0aa2)/theta2 + (sint*aa0aa2)/theta3 + (2*cost*aa[1]*aa2aa2)/theta4 + (sint*aa[1]*aa2aa2)/theta3;
        jacob[26] = (2*aa[2])/theta2 - cost*((2*aa[2])/theta2 - (2*aa2aa2aa2)/theta4) - (2*aa2aa2aa2)/theta4 + (sint*aa[2]*(aa2aa2/theta2 - T(1)))/theta;
    }
}


template<typename T, template<typename, int> class _Storage>
template<typename Other>
AngleAxis<T, _Storage>& AngleAxis<T, _Storage>::operator =(const Other& other)
{
    T m[3];
    other.template cast<T>().toAngleAxis(m);
    memcpy(Storage::data(), m, Traits<Self>::kDataSize * sizeof(T));
    return *this;
}

template<typename T, template<typename, int> class _Storage>
typename AngleAxis<T, _Storage>::OwnedSelf AngleAxis<T, _Storage>::inverse() const
{
    AngleAxis<T, OwnedStorage> inversed(*this);
    inversed.data()[0] = -inversed.data()[0];
    inversed.data()[1] = -inversed.data()[1];
    inversed.data()[2] = -inversed.data()[2];
    return inversed;
}

template<typename T, template<typename, int> class _Storage>
typename AngleAxis<T, _Storage>::OwnedSelf AngleAxis<T, _Storage>::Identity()
{
    OwnedSelf ret;
    ret.data()[0] = T(0);
    ret.data()[1] = T(0);
    ret.data()[2] = T(0);
    return ret;
}
template<typename T, template<typename, int> class _Storage>
void AngleAxis<T, _Storage>::toAngleAxis(T* angleAxis, T* jacob) const
{
    memcpy(angleAxis, Storage::data(), sizeof(T) * Storage::dataSize());
    if (!jacob) return;
    for (int i = 0; i < 3 * 3; ++i)
        jacob[i] = (i % 4 == 0) ? T(1) : T(0);
}

template<typename Derived>
PoseBase<Derived>::PoseBase()
: Storage()
, orientation_(Storage::data())
, location_(Storage::data() + kOrientSize)
{}

template<typename Derived>
PoseBase<Derived>::PoseBase(T* data)
: Storage(data)
, orientation_(Storage::data())
, location_(Storage::data() + kOrientSize)
{}

template<typename Derived>
PoseBase<Derived>::PoseBase(const T* data)
: Storage(data)
, orientation_(Storage::data())
, location_(Storage::data() + kOrientSize)
{}

template<typename T,
    template<typename, template<typename, int> class> class _Orientation,
    template<typename, int> class _Storage>
Pose<T, _Orientation, _Storage>::Pose(const Orientation& orient, const T* loc)
: BaseType()
{
    CheckOwnedStorage<_Storage> this_constructor_is_only_for_OwnedStoraged_based_object;
    BaseType::orientation_ = orient;
    BaseType::location_[0] = loc[0];
    BaseType::location_[1] = loc[1];
    BaseType::location_[2] = loc[2];
}

template<typename Derived>
void PoseBase<Derived>::toMatrix3x4(T* matrix, T* jacob) const
{
    if (!jacob)
        orientation().toRotationMatrix(matrix);
    else
    {
        orientation().toRotationMatrix(matrix, jacob);
    }
    matrix[9] = location()[0];
    matrix[10] = location()[1];
    matrix[11] = location()[2];
}

template<typename Derived>
typename PoseBase<Derived>::OwnedDerived PoseBase<Derived>::inverse() const
{
    OwnedDerived inversed;
    inversed.orientation() = orientation().inverse();
    orientation().transformIn(location(), inversed.location());
    inversed.location()[0] *= T(-1);
    inversed.location()[1] *= T(-1);
    inversed.location()[2] *= T(-1);

    return inversed;
}

template<typename Derived>
typename PoseBase<Derived>::OwnedDerived PoseBase<Derived>::Identity()
{
    OwnedDerived ret;
    ret.orientation() = Traits<Orientation>::OwnedType::Identity();
    ret.location()[0] = T(0);
    ret.location()[1] = T(0);
    ret.location()[2] = T(0);
    return ret;
}
template<typename Derived>
typename PoseBase<Derived>::OwnedDerived PoseBase<Derived>::combine(const Derived& other) const
{
    Derived combined;
    combined.orientation() = orientation().combine(other.orientation());
    orientation().transformOut(other.location(), combined.location());
    combined.location()[0] += location()[0];
    combined.location()[1] += location()[1];
    combined.location()[2] += location()[2];
    return combined;
}

template<typename T,
    template<typename, template<typename, int> class> class _Orientation,
    template<typename, int> class _Storage>
template<typename Other>
Pose<T, _Orientation, _Storage>& Pose<T, _Orientation, _Storage>::operator =(const Other& other)
{
    BaseType::orientation() = other.orientation();
    BaseType::location()[0] = other.location()[0];
    BaseType::location()[1] = other.location()[1];
    BaseType::location()[2] = other.location()[2];
    return *this;
}

template<typename T, template<typename, int> class Storage>
std::ostream& operator <<(std::ostream& os, const RotationMatrix<T, Storage>& R)
{
    os << "\t" << R.data()[0] << "\t" << R.data()[3] << "\t" << R.data()[6] << std::endl;;
    os << "\t" << R.data()[1] << "\t" << R.data()[4] << "\t" << R.data()[7] << std::endl;;
    os << "\t" << R.data()[2] << "\t" << R.data()[5] << "\t" << R.data()[8] << std::endl;;
    return os;
}

template<typename T, template<typename, int> class Storage>
std::ostream& operator <<(std::ostream& os, const AngleAxis<T, Storage>& angleAxis)
{
    os << "\t" << angleAxis.data()[0] << "\t" << angleAxis.data()[1] << "\t" << angleAxis.data()[2] << std::endl;;
    return os;
}

template<typename T,
    template<typename, template<typename, int> class> class _Orientation,
    template<typename, int> class _Storage>
std::ostream& operator <<(std::ostream& os, const Pose<T, _Orientation, _Storage>& P)
{
    os << P.orientation();
    os << "\t" << P.location()[0] << "\t" << P.location()[1] << "\t" << P.location()[2] << std::endl;
    return os;
}

typedef RotationMatrix<double, OwnedStorage> RotationMatrixd;
typedef AngleAxis<double, OwnedStorage> AngleAxisd;
typedef Pose<double, AngleAxis, OwnedStorage> Pose6d;
typedef Pose<double, RotationMatrix, OwnedStorage> Pose34d;

typedef RotationMatrix<double, SharedStorage> AsRotationMatrixd;
typedef AngleAxis<double, SharedStorage> AsAngleAxisd;
typedef Pose<double, AngleAxis, SharedStorage> AsPose6d;
typedef Pose<double, RotationMatrix, SharedStorage> AsPose34d;

typedef RotationMatrix<float, OwnedStorage> RotationMatrixf;
typedef AngleAxis<float, OwnedStorage> AngleAxisf;
typedef Pose<float, AngleAxis, OwnedStorage> Pose6f;
typedef Pose<float, RotationMatrix, OwnedStorage> Pose34f;

typedef RotationMatrix<float, SharedStorage> AsRotationMatrixf;
typedef AngleAxis<float, SharedStorage> AsAngleAxisf;
typedef Pose<float, AngleAxis, SharedStorage> AsPose6f;
typedef Pose<float, RotationMatrix, SharedStorage> AsPose34f;

} // namespace geometry 
} // namespace idl

#endif // POSE_H_
