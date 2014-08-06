#ifndef __itkConvexHull_h
#define __itkConvexHull_h

#include "itkHalfConvexHullTreeContainer.h"

namespace itk
{
template< class TVector >
class ConvexHull : public Object
{
/** \class ConvexHull
 * \brief Dynamic convex hull structure of a set of points.
 *
 * This class implements a dynamic convex hull of a set of points.
 * It uses the class HalfConvexHullTreeNode to maintain the hull.
 * It is split into two half hulls, and updated in \f$O(\log^2n)\f$.
 *
 * It does a 2D convex hull, but we can use n dimensional vectors,
 * telling in which dimensions to work.
 *
 */
public:
	/** Standard typedefs */
	typedef ConvexHull									Self;
	typedef SmartPointer< Self >						Pointer;
	typedef SmartPointer< const Self >					ConstPointer;
	typedef Object										Superclass;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( ConvexHull, Object );

	/** Some convenient tyepdefs. */
	typedef	TVector									VectorType;
	typedef typename VectorType::ValueType			ValueType;
	typedef HalfConvexHullTreeNode< VectorType >	HalfConvexHullTreeNodeType;
	typedef HalfConvexHullTreeContainer< HalfConvexHullTreeNodeType > HalfConvexHullTreeContainerType;
	typedef typename HalfConvexHullTreeNodeType::ConcatenableQueueType ConcatanableQueueType;
	typedef typename HalfConvexHullTreeNodeType::ConcatenableQueueNodeType ConcatenableQueueNodeType;

	typedef typename ConcatenableQueueNodeType::ChildIdentifier	ChildIdentifier;

	/** Left and Right indexes. */
	static const ChildIdentifier Right = ConcatenableQueueNodeType::Right;
	static const ChildIdentifier Left = ConcatenableQueueNodeType::Left;

	/** Dimension of the vector used. */
	static const unsigned int Dimension = VectorType::Dimension;

	/** Set/Get macros for the sort dimension in the binary trees. */
	bool SetDimension( const unsigned int dimension );
	itkGetConstMacro( Dimension, unsigned int );

	/** Set/Get the other dimension of the 2D plane. */
	bool SetOtherDimension( const unsigned int dimension );
	itkGetConstMacro( OtherDimension, unsigned int );

	/** Get the upper hull. */
	HalfConvexHullTreeContainerType *GetUpperHull() const { return m_UpperHull; }

	/** Get the lower hull. */
	HalfConvexHullTreeContainerType *GetLowerHull() const { return m_LowerHull; }

	/** Get macro to get the vertical width. */
	itkGetConstMacro( VerticalWidth, double );

	/** Add a point in \f$O(log^2n)\f$. */
	void Add( VectorType & vector );

	/** Remove a point of the set of points in \f$O(log^2n)\f$. */
	void Remove( VectorType & vector );

	/** Update the vertical width of the convex hull in, saying on which
	 * point the maximum vertical width is. 
	 * See "An elementary algorithm for digital line recognition 
	 * in the general case", L. Buzer, 2005.
	 */
	void UpdateVerticalWidth( const VectorType & vector );
	
	/** Update the vertical width of the convex hull. */
	void UpdateVerticalWidth();

protected:
	ConvexHull();
	virtual ~ConvexHull(){};

	/** Dimensions definition to have a plan. */
	unsigned int m_Dimension, m_OtherDimension;

	/** Vertical width of the convex hull. */
	double m_VerticalWidth;

	/** Pointers to the upper and lower hull. */
	typename HalfConvexHullTreeContainerType::Pointer m_UpperHull, m_LowerHull;

private:
	ConvexHull( const Self & );	// purposely not implemented
	void operator=( const Self & );			// purposely not implemented


};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConvexHull.hxx"
#endif

#endif
