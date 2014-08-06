#ifndef __itkLineSegmentation_h
#define __itkLineSegmentation_h

#include <itkVector.h>
#include <itkIndex.h>

namespace itk
{

template< class TIndexArray >
class LineSegmentation
{
public:
	typedef TIndexArray 										IndexArrayType;
	typedef IndexArrayType*										IndexArrayPointer;
	typedef const IndexArrayType*								IndexArrayConstPointer;
	typedef typename IndexArrayType::value_type 				IndexType;
	typedef typename IndexType::OffsetType						OffsetType;
	typedef Index< 2 > 											Index2DType;
	typedef typename Index2DType::OffsetType					Offset2DType;
	typedef int 												ValueType;

	static const unsigned int Dimension = IndexType::Dimension;
	static const unsigned int Projection2DNumber = Dimension * ( Dimension - 1 ) / 2;

	typedef Vector< ValueType, Dimension > 						VectorType;
	typedef Vector< ValueType, 2 >								Vector2DType;
	typedef FixedArray< IndexType, Dimension >					IndexFixedArrayType;

	typedef typename IndexArrayType::iterator					IndexIterator;
	typedef typename IndexArrayType::const_iterator				IndexConstIterator;

	typedef struct
	{
		VectorType MainVector;
		VectorType Offset;
		VectorType Thickness;
		unsigned int MainDimension;
		IndexConstIterator Begin, End;
		double Width;
		unsigned int Size;
	} 															DiscreteLineType;
	typedef struct
	{
		Offset2DType MainVector;
		ValueType Offset;
		ValueType Thickness;	
		unsigned int Octant;
		IndexType L, U;

		unsigned int X, Y;
		bool Lower;
		unsigned int Size;
		double Width;
	}															ProjectedBlurredNaiveLine2DType;

	typedef FixedArray< DiscreteLine2DType, ProjectionNumber >	ProjectedDiscreteLine2DFixedArrayType;
	typedef FixedArray< unsigned int, Projection2DNumber >		Projection2DFixedArrayType;

	typedef typename std::vector< DiscreteLineType >		 	DiscreteLineArrayType;
	typedef DiscreteLineArrayType*								DiscreteLineArrayPointer;
	typedef const DiscreteLineArrayType*						DiscreteLineArrayConstPointer;
	
	typedef typename DiscreteLineArrayType::iterator			DiscreteLineIterator;
	typedef typename DiscreteLineArrayType::const_iterator		DiscreteLineConstIterator;

	typedef typename std::vector< IndexConstIterator >			IndexConstIteratorContainerType;
	typedef typename std::vector< IndexIterator >				IndexIteratorContainerType;

	DiscreteLineArrayPointer GetSegmentedLine() { return m_SegmentedLine; }
	void SetInput( IndexArrayConstPointer input ) { m_Input = input; }
	IndexArrayConstPointer GetInput() { return m_Input; }

	virtual DiscreteLineArrayPointer SegmentLine();

	LineSegmentation(){ m_SegmentedLine = new DiscreteLineArrayType; }
	~LineSegmentation(){};

protected:
	

	inline void SelectSegmentedLines();
	inline void Segment( IndexConstIterator begin, IndexConstIterator end);
	inline bool UpdateLU( IndexConstIteratorContainerType & LF, IndexConstIteratorContainerType & UF, IndexConstIteratorContainerType & LL, IndexConstIteratorContainerType & UL, const DiscreteLineType & line, IndexConstIterator it, IndexConstIterator jt ) const;
	virtual DiscreteLineType ComputeDiscreteLine( IndexConstIteratorContainerType & LF, IndexConstIteratorContainerType & UF, IndexConstIteratorContainerType & LL, IndexConstIteratorContainerType & UL, IndexConstIterator begin, IndexConstIterator it, const bool & lower, unsigned int & dimension ) const;

	inline bool IsOnBoundary( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const;
	inline bool IsOnBoundary( const IndexType & M, const DiscreteLineType & line ) const;
	inline bool IsWeaklyExterior( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const;
	inline bool IsWeaklyExterior( const IndexType & M, const DiscreteLineType & line ) const;
	inline bool IsExterior( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const;
	inline bool IsExterior( const IndexType & M, const DiscreteLineType & line ) const;
	inline bool IsStronglyExterior( const IndexType & M, const DiscreteLineType & line, bool & lower, unsigned int & dimension ) const;
	inline bool IsStronglyExterior( const IndexType & M, const DiscreteLineType & line ) const;

	inline ValueType Cost( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline ValueType Cost( const Index2DType & idx, const Vector2DType & vector ) const;
	inline ValueType Cost( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;

	inline bool IsIncreasing( const DiscreteLineType & line, unsigned int & i ) const;
	inline bool IsDecreasing( const DiscreteLineType & line, unsigned int & i ) const;

	inline bool IsExteriorLower( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsExteriorLower( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsExteriorUpper( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsExteriorUpper( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsExterior( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsExterior( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsWeaklyExteriorLower( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsWeaklyExteriorLower( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsWeaklyExteriorUpper( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsWeaklyExteriorUpper( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsWeaklyExterior( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsWeaklyExterior( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsOnLowerBoundary( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsOnLowerBoundary( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsOnUpperBoundary( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsOnUpperBoundary( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsOnBoundary( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsOnBoundary( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsStronglyExteriorLower( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsStronglyExteriorLower( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsStronglyExteriorUpper( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsStronglyExteriorUpper( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsStronglyExterior( const ValueType & value, const DiscreteLine2DType & line ) const;
	inline bool IsStronglyExterior( const Index2DType & idx, const DiscreteLine2DType & line ) const;
	inline bool IsExteriorLower( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsExteriorLower( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsExteriorUpper( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsExteriorUpper( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsExterior( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsExterior( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsWeaklyExteriorLower( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsWeaklyExteriorLower( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsWeaklyExteriorUpper( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsWeaklyExteriorUpper( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsWeaklyExterior( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsWeaklyExterior( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsOnLowerBoundary( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsOnLowerBoundary( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsOnUpperBoundary( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsOnUpperBoundary( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsOnBoundary( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsOnBoundary( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsStronglyExteriorLower( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsStronglyExteriorLower( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsStronglyExteriorUpper( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsStronglyExteriorUpper( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsStronglyExterior( const ValueType & value, const DiscreteLineType & line, unsigned int i ) const;
	inline bool IsStronglyExterior( const IndexType & idx, const DiscreteLineType & line, unsigned int i ) const;


private:

	IndexArrayConstPointer m_Input;
	DiscreteLineArrayPointer m_SegmentedLine;

};

}// namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLineSegmentation.hxx"
#endif

#endif
