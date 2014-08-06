#ifndef __itkFiberExtractionGraphFilter_h
#define __itkFiberExtractionGraphFilter_h

#include "Graph/itkGraphToGraphFilter.h"

namespace itk
{
/** \class FiberExtractionGraphFilter
 * \brief Graph filter that transforms an input skeleton graph by merging
 * the nodes of it depending of the connections between lines.
 *
 * Depending of the given TGraph, that has to be herited from SkeletonGraph,
 * the output graph can be different. This class is just a mechanism that
 * calls functions from the graph, that heps decide which nodes should be
 * merged or not.
 *
 */
template< class TGraph >
class ITK_EXPORT FiberExtractionGraphFilter : public GraphToGraphFilter< TGraph, TGraph >
{
public:
	/** Standard class typedefs. */
	typedef FiberExtractionGraphFilter				Self;
	typedef GraphToGraphFilter< TGraph, TGraph >	Superclass;
	typedef SmartPointer< Self >					Pointer;
	typedef SmartPointer< const Self >				ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( GraphToGraphFilter, FiberExractionGraphFilter );

	/** Some convenient typedefs. */
	typedef TGraph													GraphType;
	typedef typename GraphType::Pointer								GraphPointer;
	typedef typename GraphType::GraphTraitsType						GraphTraitsType;
	typedef typename GraphType::NodeContainerType					NodeContainerType;
	typedef typename GraphType::EdgeContainerType					EdgeContainerType;
	typedef typename GraphTraitsType::NodeType						NodeType;
	typedef typename GraphTraitsType::NodePointerType				NodePointerType;
	typedef typename GraphTraitsType::EdgeType						EdgeType;
	typedef typename GraphTraitsType::EdgePointerType				EdgePointerType;
	typedef typename GraphTraitsType::VectorType					VectorType;
	typedef typename GraphTraitsType::IndexContainerType			IndexContainerType;
	typedef typename GraphTraitsType::NodeIdentifierType			NodeIdentifierType;
	typedef typename GraphTraitsType::EdgeIdentifierType			EdgeIdentifierType;
	typedef typename IndexContainerType::value_type					IndexType;
	typedef typename GraphTraitsType::ScalarContainerType			ScalarContainerType;
	typedef typename GraphTraitsType::NodeWeightType				NodeWeightType;
	typedef typename GraphTraitsType::EdgeIdentifierContainerType	EdgeIdentifierContainerType;

	/** Definition of the dimension of the represented lines. */
	static const unsigned int Dimension = IndexType::Dimension;


	/** Set/Get the threshold that allow two nodes to be merged. */
	itkGetConstMacro( MergeThreshold, double );
	itkSetMacro( MergeThreshold, double );

	/** Set/Get the split that says when to split a line in 2 parts. */
	itkGetConstMacro( SplitThreshold, double );
	itkSetMacro( SplitThreshold, double );

	/** Set/Get max intersection identifier of the graph (FrontLabel and BackLabel). */
	itkGetConstMacro( MaxIntersection, NodeIdentifierType );
	void SetMaxIntersection( NodeIdentifierType inter );

protected:

	FiberExtractionGraphFilter();
	~FiberExtractionGraphFilter(){};

	/** GenerateData function. */
	virtual void GenerateData();

	/** GenerateOutputInformation function. */
	virtual void GenerateOutputInformation();
	
	/** Initialize the input graph. */
	virtual void Initialize();

	/** Split and merge routine. */
	virtual void SplitAndMerge();

	NodeIdentifierType m_MaxIntersection;
	double m_MergeThreshold, m_SplitThreshold;
	NodeWeightType m_Maximum;
	bool m_MaxIntersectionSet;

private:
	FiberExtractionGraphFilter( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely node implemented

	
	GraphPointer m_Output;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFiberExtractionGraphFilter.hxx"
#endif

#endif
