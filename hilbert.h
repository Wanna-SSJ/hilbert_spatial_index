#ifndef HILBERT_H_
#define HILBERT_H_
#include <cstdint>
#include <istream>
#include <stdexcept>
class GeoDegreeCoordinateSystem {
public:
  static void GetRanges(double xrange[2], double yrange[2]) {
    xrange[0] = -180;
    xrange[1] = 180;
    yrange[0] = -90;
    yrange[1] = 90;
  }
};

template <typename CoordinateSystem>
class GeoPoint {
    private:
        double x_;
        double y_;

    public:
        GeoPoint():x_(0),y_(0) {}
        GeoPoint(double x, double y) : x_(x), y_(y) {}

        GeoPoint(const GeoPoint<CoordinateSystem>& point) : x_(point.x_), y_(point.y_) {}

        GeoPoint<CoordinateSystem>& operator=(const GeoPoint<CoordinateSystem>& point) {
            x_ = point.x_;
            y_ = point.y_;
            return *this;
        }

        bool operator==(const GeoPoint<CoordinateSystem>& point) const {
            return x_ == point.x_ && y_ == point.y_;
        }

        // Get the x coordinate.
        double x() const {
            return x_;
        }

        // Get the y coordinate.
        double y() const {
            return y_;
        }

        // For pretty printing.
        friend std::ostream& operator<<(std::ostream& out, const GeoPoint<CoordinateSystem>& point) {
            out << "(" << point.x_ << "," << point.y_ << ")";
            return out;
        }
};

// GeoHash wraps a HashType geohash value within CoordinateSystem.
template <typename CoordinateSystem, typename HashType>
class GeoHash {
    private:
      HashType value_;
    public:
      GeoHash(HashType value) : value_(value) {}

      GeoHash(const GeoHash& hash) : value_(hash.value_) {}

      GeoHash<CoordinateSystem, HashType>& operator=(const GeoHash<CoordinateSystem, HashType>& hash) {
        value_ = hash.value_;
        return *this;
      }

      bool operator==(const GeoHash<CoordinateSystem, HashType>& hash) const {
        return value_ == hash.value_;
      }

      // Get the hash value.
      HashType value() const {
        return value_;
      }
};

// Used for mapping from hash type to coordinate type and for compile-time validation of the hash
// type:
// A 64 bit hash type carries 32 bits from each coordinate.
// A 32 bit hash type carries 16 bits from each coordinate.
// A 16 bit hash type carries 8 bits from each coordinate.
template <typename HashType>
struct HashTypeInterpreter;

template <>
struct HashTypeInterpreter<uint64_t> {
  typedef uint32_t CoordinateType;
};

template <>
struct HashTypeInterpreter<uint32_t> {
  typedef uint16_t CoordinateType;
};

template <>
struct HashTypeInterpreter<uint16_t> {
  typedef uint8_t CoordinateType;
};

// Maps and unmaps points to distances down the Hilbert curve of #bits(sizeof(CoordinateType))
// iterations.
//
// To map a point to its Hilbert value, this implementation operates iteratively by dividing the
// current area into quadrants, calculating the quadrant into which the point falls, and rotating
// the quadrant ordering to match the current Hilbert curve iteration. This process repeats with the
// calculated quadrant for #bits(sizeof(CoordinateType)) iterations.
//
// To map a Hilbert value to a point, this implementation operates iteratively by extracting the
// rotated quadrant from the appropriate two-bits of the Hilbert value and un-rotating based upon
// the Hilbert curve iteration which produced the bits. This process repeats for
// #bits(sizeof(CoordinateType)) iterations.
//
// See http://en.wikipedia.org/wiki/Hilbert_curve for the construction of a Hilbert curve.
template <typename HashType>
class Hilbert {
private:
  typedef typename HashTypeInterpreter<HashType>::CoordinateType CoordinateType;

public:
  enum Rotation {Up, Right, Down, Left};

  HashType hilbert(CoordinateType x, CoordinateType y) const {
    return hilbert(x, y, nullptr);
  }

  // Map from the point (x, y) to the length down the Hilbert curve of #bits(sizeof(CoordinateType))
  // iterations.
  // rotations: #bits(sizeof(CoordinateType)) Rotation wide array which will be populated with the
  //            quadrant rotations at each iteration of hilbert curve generation; ignored if
  //            nullptr.
  HashType hilbert(CoordinateType x, CoordinateType y, Rotation* rotations) const {
    HashType result = 0;
    Rotation rotation = Up;
    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      uint8_t quadrant = 0;
      // | 0 | 3 |
      //  -------
      // | 1 | 2 |
      if (x & (1 << i)) {
        if (y & (1 << i)) {
          quadrant = 3;
        } else {
          quadrant = 2;
        }
      } else {
        if (y & (1 << i)) {
          quadrant = 0;
        } else {
          quadrant = 1;
        }
      }

      if (rotations != nullptr) {
        rotations[i] = rotation;
      }
      HashType rotated = rotate(rotation, quadrant);
      rotation = new_rotation(rotation, rotated);
      result |= (rotated << (i << 1));
    }
    return result;
  }

  void unhilbert(HashType d, CoordinateType& x, CoordinateType& y) const {
    unhilbert(d, x, y, nullptr);
  }

  // Map from d, the output of the hilbert function, to the point (x, y).
  // rotations: #bits(sizeof(CoordinateType)) Rotation wide array which will be populated with the
  //            quadrant rotations at each iteration of hilbert curve generation; ignored if
  //            nullptr.
  void unhilbert(HashType d, CoordinateType& x, CoordinateType& y, Rotation* rotations) const {
    x = 0;
    y = 0;
    Rotation rotation = Up;
    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      if (rotations != nullptr) {
        rotations[i] = rotation;
      }
      uint8_t quadrant = ((d & (3ll << (i << 1))) >> (i << 1));
      uint8_t unrotated = rotate(rotation, quadrant);
      rotation = new_rotation(rotation, quadrant);
      // | 0 | 3 |
      //  -------
      // | 1 | 2 |
      switch (unrotated) {
        case 0:
          y |= (1 << i);
          break;
        case 1:
          // no op
          break;
        case 2:
          x |= (1 << i);
          break;
        case 3:
          x |= (1 << i);
          y |= (1 << i);
          break;
        default:
          throw std::runtime_error("internal error");
      }
    }
  }

private:
  // Rotate quadrant by rotation. The provided quadrant is with respect to the Up rotation.
  //
  // Note that this function is self-inverting in the sense that rotate(r, rotate(r, x)) == x.
  uint8_t rotate(Rotation rotation, uint8_t quadrant) const {
    // | 2 | 3 |
    //  -------
    // | 1 | 0 |
    static uint8_t right[4] = {2, 1, 0, 3};
    // | 2 | 1 |
    //  -------
    // | 3 | 0 |
    static uint8_t down[4] = {2, 3, 0, 1};
    // | 0 | 1 |
    //  -------
    // | 3 | 2 |
    static uint8_t left[4] = {0, 3, 2, 1};

    switch (rotation) {
      case Up:
        return quadrant;
      case Right:
        return right[quadrant];
      case Down:
        return down[quadrant];
      case Left:
        return left[quadrant];
      default:
        throw std::runtime_error("internal error");
    }
  }

  // Compute the new rotation for quadrant and rotation. The provided quadrant is with respect to
  // rotation.
  Rotation new_rotation(Rotation rotation, uint8_t quadrant) const {
    Rotation result = rotation;
    // | Left | Right |
    //  --------------
    // |  Up  |   Up  |
    if (quadrant == 0) {
      switch (rotation) {
        case Up:
          result = Left;
          break;
        case Right:
          result = Down;
          break;
        case Down:
          result = Right;
          break;
        case Left:
          result = Up;
          break;
      }
    } else if (quadrant == 3) {
      switch (rotation) {
        case Up:
          result = Right;
          break;
        case Right:
          result = Up;
          break;
        case Down:
          result = Left;
          break;
        case Left:
          result = Down;
          break;
      }
    }
    return result;
  }
};

// GeoHasher hashes and unhashes GeoPoints. The GeoHash produced by geohashing a point is computed
// as the distance down the Hilbert curve of #bits(sizeof(CoordinateType)) iterations overlaid onto
// the CoordinateSystem's range. Each coordinate is first mapped to a #bits(sizeof(CoordinateType))
// bit integer which represents the coordinate's location within the CoordinateSystem's range as
// follows: The high-order bit gives the coordinate's location relative to the midpoint of the
// range, the next bit gives the coordinate's location relative to the the midpoint of the relevant
// half of the range, and so on. The results of this process are the coordinates within an integer
// coordinate system that spans the CoordinateSystem's range onto which a Hilbert curve can be
// overlaid.
//
// The construction of the Hilbert curve gives the nice property that points which are near one
// another spatially tend to fall near one another on the curve; hence, points that are near one
// another spatially tend to produce GeoHashes that are close in magnitude.
template <typename HashType>
class GeoHasher {
private:
  typedef typename HashTypeInterpreter<HashType>::CoordinateType CoordinateType;
private:
  Hilbert<HashType> hilbert_;
public:
  // Compute the GeoHash for point.
  template <typename CoordinateSystem>
  GeoHash<CoordinateSystem, HashType> hash(const GeoPoint<CoordinateSystem>& point) const {
    double x_range[2];
    double y_range[2];
    CoordinateSystem::GetRanges(x_range, y_range);

    CoordinateType x = 0;
    CoordinateType y = 0;
    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      double y_mid = (y_range[0] + y_range[1]) / 2.;
      //std::cout << x << ',' << y <<std::endl;
      if (point.y() >= y_mid) {
        y |= 1 << i;
        y_range[0] = y_mid;
      } else {
        y_range[1] = y_mid;
      }
      double x_mid = (x_range[0] + x_range[1]) / 2.;
      if (point.x() >= x_mid) {
        x |= 1 << i;
        x_range[0] = x_mid;
      } else {
        x_range[1] = x_mid;
      }
    }
    return hilbert_.hilbert(x, y);
  }

  // Compute the GeoPoint for hash.
  template <typename CoordinateSystem>
  GeoPoint<CoordinateSystem> unhash(const GeoHash<CoordinateSystem, HashType>& hash) const {
    double x_range[2];
    double y_range[2];
    CoordinateSystem::GetRanges(x_range, y_range);

    CoordinateType x = 0;
    CoordinateType y = 0;
    hilbert_.unhilbert(hash.value(), x, y);

    for (int i = (sizeof(CoordinateType) << 3) - 1; i >= 0; --i) {
      double y_mid = (y_range[0] + y_range[1]) / 2.;
      if (y & (1 << i)) {
        y_range[0] = y_mid;
      } else {
        y_range[1] = y_mid;
      }
      double x_mid = (x_range[0] + x_range[1]) / 2.;
      if (x & (1 << i)) {
        x_range[0] = x_mid;
      } else {
        x_range[1] = x_mid;
      }
    }
    return GeoPoint<CoordinateSystem>(x_range[0], y_range[0]);
  }
};
#endif // HILBERT_H_
