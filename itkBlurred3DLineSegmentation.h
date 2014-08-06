#ifndef __itkBlurred3DLineSegmentation_h
#define __itkBlurred3DLineSegmentation_h

#include <itkVector.h>
#include <itkIndex.h>
#include "itkConvexHull.h"

namespace itk
{
/** \class Blurred3DLineSegmentation
 * \brief Tool to segment a 3D discrete line into multiple
 * 3D blurred segments.
 *
 * This class uses the concept of blurred 3D line and uses the
 * results of Isabelle Debled - Rennesson and Thanh Phuong Nguyen
 * in the paper "On the local properties of digital curves.
 *
 * See at the following link :
 * http://hal.archives-ouvertes.fr/docs/00/43/73/04/PDF/JournalShapeModeling_NguyenDebled.pdf
 */
template< class TIndexArray >
class Blurred3DLineSegmentation
{
public:

	Blurred3DLineSegmentation(){}
	~Blurred3DLineSegmentation(){};


	/** Dimension of the vectors for the blurred segment. */
	static const unsigned int Dimension = 3;

	/** Number of possible projections into 2D plans. */
	static const unsigned int ProjectionCount = 3;

	typedef Vector< double, Dimension >			VectorType;

	typedef ConvexHull< VectorType > 			ConvexHullType;
	typedef typename ConvexHullType::ConcatenableQueueNodeType ConcatenableQueueNodeType;

	/** 2D line structure. */
	typedef struct
	{
		unsigned int Octant;
		unsigned int a, b;
		int Offset;
		unsigned int Thickness;
		ConcatenableQueueNodeType *LL, *UL;
	}											Discrete2DLineType;

	typedef std::vector< Discrete2DLineType >	DiscreteLineType;


	typedef std::vector< DiscreteLineType > 	MaxDiscreteSegmentContainerType;



	/** Routine to segment the set of point from begin to end into a list of
	 * blurred segments of width m_Width. 
	 */
	template< class TIterator >
	MaxDiscreteSegmentContainerType Segment( TIterator begin, TIterator end);

	/** Update the line after an addition or a removal of a point in the set.
	 * M is the removed or added point.
	 */
	void UpdateLine( const ConvexHullType *hull, Discrete2DLineType & line, ConcatenableQueueNodeType *M, bool increase ) const;

	/** Tells if the line is a blurred segment of width m_Width or not. */
	bool IsBlurredLine( const Discrete2DLineType & line ) const;

	/** Add a point to the set of points. */
	void Add( ConvexHullType *hull, VectorType point ) const;

	/** Remoe a point to the set of points. */
	void Remove( ConvexHullType *hull, VectorType point ) const;

	/** Set/Get the width of the blurred segments. */
	double GetWidth() const { return m_Width; }
	void SetWidth( double & width ) const { m_Width = width; };

protected:

	double m_Width;


private:

};

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlurred3DLineSegmentation.hxx"
#endif

#endif
