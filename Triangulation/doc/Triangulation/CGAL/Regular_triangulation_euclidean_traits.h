
namespace CGAL {

/*!
\ingroup PkgTriangulationsTraitsClasses

The class `Regular_triangulation_euclidean_traits` is used internally by the
class `Regular_triangulation` to wrap its first template parameter
(`RegularTriangulationTraits_`)
so that the base class `Triangulation` manipulates weighted points instead 
of bare points.

\tparam K must be a model of the `RegularTriangulationTraits` concept.

In addition to the types described below, the following predicates and functors
are adapted so that they can be called
with weighted points as parameters. In practice, the functors from the
base class `K` are called, ignoring the weights.
- `Orientation_d`
- `Construct_flat_orientation_d`
- `In_flat_orientation_d`
- `Contained_in_affine_hull_d`
- `Compare_lexicographically_d`
- `Compute_coordinate_d`
- `Point_dimension_d`
- `Less_coordinate_d`

*/

template <typename K>
class Regular_triangulation_euclidean_traits : public K {
public:

  /// \name Types
  /// @{

  /*!
  The weighted point type.
  */
  typedef typename K::Weighted_point_d              Point_d;

  /*!
  The (un-weighted) point type.
  */
  typedef typename K::Point_d                       Bare_point_d;

  /// @}

  /// \name Creation
  /// @{

  /*!
  The default constructor.
  */
  Regular_triangulation_euclidean_traits();

  /// @}

};

} /* end namespace CGAL */
