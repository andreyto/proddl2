//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//#
//#   See COPYING file distributed along with the PRODDL package for the
//#   copyright and license terms.
//#
//### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
//
//
#ifndef PRODDL_TRANSFORMATION_H__
#define PRODDL_TRANSFORMATION_H__

#include "PRODDL/Math/matr.hpp"

#include "PRODDL/Geom/euler.hpp"

#include <iostream> // for debugging printouts

/*
    Linear coordinate transformation.

    This is a C++ reimplementation of ScientificPython Transformation module.

    Transformation objects represent linear coordinate transformations
    in a 3D space. They can be applied to vectors, returning another vector.
    If 't' is a transformation and 'v' is a vector, 't(v)' returns
    the transformed vector.

    Transformations support composition: if 't1' and 't2' are transformation
    objects, 't1*t2' is another transformation object which corresponds
    to applying t1 *AFTER* t2.

*/

namespace PRODDL { namespace Geom { 


  // Forward declarations

  template <typename T_num> class Translation;
  template <typename T_num> class Rotation;
  template <typename T_num> class RotationTranslation;


//
// Pure translation
//

template<typename T_num> class Translation {

public:

  typedef typename SpaceTraits<T_num>::Point3 Point;
  typedef typename SpaceTraits<T_num>::Matrix3x3 Matrix3;

  typedef typename SpaceTraits<T_num>::VPoint3 Points;

protected:

  Point vector;

public:

  const Point& getVector() const {
    return vector;
  }


  explicit Translation(const Point& _vector):
    vector(_vector)
  {}

  explicit Translation(T_num x):
    vector(x,x,x)
  {}

  Translation():
    vector(0,0,0)
  {}

  Translation operator* (const Translation& other) const {
    return Translation(vector + other.getVector());
  }

  RotationTranslation<T_num> operator* (const Rotation<T_num>& other) const {
    return RotationTranslation<T_num>(other.getTensor(),vector);
  }

  RotationTranslation<T_num> operator* (const RotationTranslation<T_num>& other) const {
    return RotationTranslation<T_num>(other.getTensor(),other.getVector()+vector);
  }

  Point operator() (const Point& _vector) const {
    return vector + _vector;
  }
  
  void operator() (Points& points) const {
    for(int i = 0; i < points.size(); i++) {
      Point& other = points(i);
      other = operator()(other);
    }
  }
  

  Point displacement() const {
    return vector;
  }


   Rotation<T_num> rotation() const {
     Matrix3 m;
     m = 
       1,0,0,
       0,1,0,
       0,0,1;
     return Rotation<T_num>(m);
   }


  Translation translation() const {
    return *this;
  }


  Translation inverse() const {
    return Translation(Point(0)-vector);
  }

};


//
// Pure rotation
//

template<typename T_num> class Rotation {

public:
  
  typedef typename SpaceTraits<T_num>::Point3 Point;
  typedef typename SpaceTraits<T_num>::Matrix3x3 Matrix3;

  typedef typename SpaceTraits<T_num>::VPoint3 Points;
  
protected:
  
  Matrix3 tensor;
  
public:
  
  
  const Matrix3& getTensor() const {
    return tensor;
  }

  explicit Rotation(const Matrix3& _tensor):
    tensor(_tensor)
  {}

  Rotation()
  {
    tensor = 
      1,0,0,
      0,1,0,
      0,0,1;
  }

  // construct from euler angles in "xyz" convention passed in Point variable
  
  explicit Rotation(const Point& euler_angles_xyz) {
    tensor = Transforms::Eu_xyz::eulerAngToRotMatr<T_num>(euler_angles_xyz);
  }
  
  Point eulerAngles() const {
    bool isUnique;
    return Transforms::Eu_xyz::rotMatrToEulerAng<T_num>(tensor,isUnique);
  }

  // bool& isUnique is the output parameter, upon return,
  // it will show if the returned triplet of Euler angles 
  // represents a unique solution.
  Point eulerAngles( bool& isUnique ) const {
    return Transforms::Eu_xyz::rotMatrToEulerAng<T_num>(tensor,isUnique);
  }

  Rotation operator* (const Rotation& other) const {
    Matrix3 m;
    m = blitz_ext::product(tensor,other.getTensor());
    return Rotation(m);
  }
  
  RotationTranslation<T_num> operator* (const Translation<T_num>& other) const {
    return RotationTranslation<T_num>(tensor,blitz_ext::product(tensor,other.getVector()));
  }
  
  RotationTranslation<T_num> operator* (const RotationTranslation<T_num>& other) const {
    //BUG:
    //Blitz does not define conversion ctor from matrixMatrixProduct expression to Matrix,
    //therefore we have to create a temporary matrix and do the assignment to it
    Matrix3 m;
    m = blitz_ext::product(tensor,other.getTensor());
    return RotationTranslation<T_num>(m,
				      blitz_ext::product(tensor,other.getVector()));
  }


  Point operator() (const Point& other) const {
    return blitz_ext::product(tensor,other);
  }
  
  Matrix3 operator() (const Matrix3& other) const {
    Matrix3 m;
    m = blitz_ext::product(other,tensor);
    Matrix3 m_ret;
    m_ret = blitz_ext::product(PRODDL::Math::inverse(tensor,1e-6),m);
    return m_ret;
  }

  void operator() (Points& points) const {
    for(int i = 0; i < points.size(); i++) {
      Point& other = points(i);
      other = operator()(other);
    }
  }

  Rotation rotation() const {
    return *this;
  }

  
  Translation<T_num> translation() const {
    return Translation<T_num>(Point(0,0,0));
  }

  Rotation inverse() const {
    return Rotation(PRODDL::Math::transpose(tensor));
  }

};

//
// Combined translation and rotation
//

template<typename T_num> class RotationTranslation {

public:
  
  typedef typename SpaceTraits<T_num>::Point3 Point;
  typedef typename SpaceTraits<T_num>::Matrix3x3 Matrix3;

  typedef typename SpaceTraits<T_num>::Point3Pair PointPair;

  typedef typename SpaceTraits<T_num>::VPoint3 Points;

protected:
  
  Matrix3 tensor;
  Point   vector;
  
public:
  
  const Matrix3& getTensor() const {
    return tensor;
  }

  const Point& getVector() const {
    return vector;
  }


  RotationTranslation(const Matrix3& _tensor, const Point& _vector) :
    tensor(_tensor),
    vector(_vector)
  {}

  RotationTranslation(const PointPair& rigidBodyCoords)
  {
    *this = Translation<T_num>(rigidBodyCoords(0)) * Rotation<T_num>(rigidBodyCoords(1));
  }

  RotationTranslation(const Point& angles, const Point& displ)
  {
    *this = Translation<T_num>(displ) * Rotation<T_num>(angles);
  }

  RotationTranslation():
    vector(0,0,0)
  {
    tensor = 
      1,0,0,
      0,1,0,
      0,0,1;
  }

  // ctors to provide conversion from Translation and Rotation objects

  RotationTranslation(const Translation<T_num>& other) {
    tensor = 
      1,0,0,
      0,1,0,
      0,0,1;
    vector = other.getVector();
  }
  
  RotationTranslation(const Rotation<T_num>& other):
    tensor(other.getTensor()),
    vector(0,0,0)
  {}

  RotationTranslation operator* (const Rotation<T_num>& other) const {
    Matrix3 m;
    m = blitz_ext::product(tensor,other.getTensor());
    return RotationTranslation<T_num>(m,
				      vector);
  }


  RotationTranslation operator* (const Translation<T_num>& other) const {
    return RotationTranslation<T_num>(tensor,
				      blitz_ext::product(tensor,other.getVector())+vector);
  }

  RotationTranslation operator* (const RotationTranslation<T_num>& other) const {
    Matrix3 m;
    m = blitz_ext::product(tensor,other.getTensor());
    Point p = blitz_ext::product(tensor,other.getVector())+vector;
    return RotationTranslation<T_num>(m,p);
  }

  Point operator() (const Point& other) const {
    return blitz_ext::product(tensor,other) + vector;
  }

  void operator() (Points& points) const {
    for(int i = 0; i < points.size(); i++) {
      Point& point = points(i);
      point = operator()(point);
    }
  }

  Rotation<T_num> rotation() const {
    return Rotation<T_num>(tensor);
  }

  
  Translation<T_num> translation() const {
    return Translation<T_num>(vector);
  }

  RotationTranslation inverse() const {
    return Rotation<T_num>(PRODDL::Math::transpose(tensor))*Translation<T_num>(0 - vector);
  }

  void anglesAndDisplacement(Point& angles,Point& displ) const {

    angles = rotation().eulerAngles();

    displ = translation().displacement();

  }

};


  // Typedefs for arrays of transformations

  template<typename T_num> struct TransformationTraits {
    typedef typename PRODDL::common_types::num_vector_type<Translation<T_num> >::Type VTranslation;
    typedef typename PRODDL::common_types::num_vector_type<Rotation<T_num> >::Type VRotation;
    typedef typename PRODDL::common_types::num_vector_type<RotationTranslation<T_num> >::Type VRotationTranslation;
  };  


  template<typename T_num>
  std::ostream& operator<< (std::ostream& out, const Translation<T_num>& x) {

    out << "Translation( " << x.getVector() << " )";
    return out;
  }

  template<typename T_num>
  std::ostream& operator<< (std::ostream& out, const Rotation<T_num>& x) {

    out << "Rotation( " << x.eulerAngles() << " )";
    return out;
  }

  template<typename T_num>
  std::ostream& operator<< (std::ostream& out, const RotationTranslation<T_num>& x) {

    out << "RotationTranslation( " << x.translation() << ", " << x.rotation() << " )";
    return out;
  }

}} // namespace PRODDL::Geom

#endif // PRODDL_TRANSFORMATION_H__
