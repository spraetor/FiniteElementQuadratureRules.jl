// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_QUADRATURE_${name | uppercase}_HH
#define DUNE_GEOMETRY_QUADRATURE_${name | uppercase}_HH

#ifndef DUNE_INCLUDING_IMPLEMENTATION
#error This is a private header that should not be included directly.
#error Use #include <dune/geometry/quadraturerules.hh> instead.
#endif

namespace Dune {

  /************************************************
   * Quadraturerule for ${basicType}
   * Generated from dune-quadrature on ${date}
   *************************************************/

  /** \brief Quadrature rules for ${basicType}
      \ingroup Quadrature
   */
  template<typename ct, int dim>
  class ${name}QuadratureRule;

  /** \brief Quadrature rules for ${basicType}
      \ingroup Quadrature
   */
  template<typename ct>
  class ${name}QuadratureRule<ct, ${dim}> : public QuadratureRule<ct, ${dim}>
  {
  public:
    /** \brief The highest quadrature order available */
    enum { highest_order = ${maxdegree} };
  private:
    friend class QuadratureRuleFactory<ct, ${dim}>;
    ${name}QuadratureRule (int p);
    ~${name}QuadratureRule() {}
  };

  template<typename ct>
  ${name}QuadratureRule<ct, ${dim}>::${name}QuadratureRule(int order)
    : QuadratureRule<ct, ${dim}>(GeometryType(GeometryType::${basicType}, ${dim}))
  {
    switch (order)
    {
    %for rule in rules:
      case ${rule['degree']}:
      %if rule and rule['coordinates'] and len(rule['coordinates']) > 0:
<%      points, weights = rule['coordinates'], rule['weights'] %>\
        // dim = ${dim}, order = ${rule['degree']}, npoints = ${len(points)}, properties = ${rule['properties']}
        % if len(rule['reference']) > 0:
        // Source: ${rule['reference']}
        % endif
        this->delivered_order = ${rule['degree']};
        this->resize(${len(points)});
        %for i,p in enumerate(points):
        (*this)[${i}] = {
          { ${(",\n" + " "*12).join(map(cast, map(str, p)))} },
          ${str(weights[i]) | cast} };
        %endfor
        break;

      %endif
    %endfor
      default:
        DUNE_THROW(QuadratureOrderOutOfRange,
                  "QuadratureRule for order " << order << " and GeometryType "
                                              << this->type() << " not available");
    }
  }

} // end namespace Dune

#endif // DUNE_GEOMETRY_QUADRATURE_${name | uppercase}_HH
