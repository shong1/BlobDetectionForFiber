#ifndef __itkSkeletonGraph_h
#define __itkSkeletonGraph_h

#include "Graph/itkGraph.h"

namespace itk
{
/** \class SkeletonGraph
 * \brief Skeleton graph data structure. Can be created from SkeletonToGraphFilter
 * 
 * The nodes of the graph represent lines, and edges represent the existence of an
 * intersection between two lines (nodes). This is a base class to represent
 * a graph of connected lines. Some functions calculating node and edge attributes
 * can be overload if attributes are different.
 *
 * This class assumes that a node (a line), is equivalent to a point, a director
 * vector, an average distance from the vector, a force (geometrical), and a
 * mean vesselness (grayscale likelihood). Herited classes can remove those, 
 * or add some, depending of TGraphTraits definition. A good template here is
 * LineGraphTraits. Edges can be represented as the concatenation of two nodes (lines).
 * Attributes are then the point, the director vector, the average distance to
 * the vector, and a float f, representing the straigness between the two lines.
 *
 * Any herited class from this base class need to have BackLabel and FrontLabel
 * attributes in the node definitin in TGraphTraits. It is used to get the
 * side of the connection between two lines.
 *
 * This kind of graph is built so we can recognize lines in a skeleton that
 * would represent a gallery of connected lines.
 * The principle here is to calculate a grayscale likelihood and a geometrical one.
 * Let call v the grayscale likelyhood and F the geometrical one.
 * \f$ v \times F \f$ is a good estimation of the likelihood to be a fiber.
 * Depending of the geometrical model, we can build a more or less
 * selective likelihood. This class considers a line as a vector and a point.
 * Thus, curved lines shall not be recognized here. Redefining a geometrical
 * likelihood might get those.
 *
 */
template< class TGraphTraits >
class SkeletonGraph : public Graph< TGraphTraits >
{
public:
	/** Standard class typedefs. */
	typedef SkeletonGraph				Self;
	typedef Graph< TGraphTraits >		Superclass;
	typedef SmartPointer< Self >		Pointer;
	typedef SmartPointer< const Self >	ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( Graph, SkeletonGraph );

	/** Some convenient typedefs. */
	typedef TGraphTraits											GraphTraitsType;
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
	typedef typename Superclass::NodeContainerType					NodeContainerType;
	typedef typename Superclass::EdgeContainerType					EdgeContainerType;
	typedef typename GraphTraitsType::NodeWeightType				NodeWeightType;
	typedef typename GraphTraitsType::EdgeWeightType				EdgeWeightType;
	typedef typename GraphTraitsType::EdgeIdentifierContainerType	EdgeIdentifierContainerType;
	typedef typename IndexType::OffsetType							OffsetType;
	typedef typename std::vector< OffsetType >						OffsetContainerType;

	/** Definition of the dimension of the represented lines. */
	static const unsigned int Dimension = IndexType::Dimension;

	/** Creates a copy of the graph. */
	virtual inline Pointer Copy() const;

	/** Initializing graph attributes. */
	virtual inline void InitializeGraphAttributes();

	/** Initializing edges from the sight of outgoing edges of a node.
	 * The node should already be initialized.
	 */
	virtual inline void InitializeOutgoingEdgesInteraction( const NodeType & node );

	/** Initializing node interaction with its neighbors (after a first self
	 * initialization of the nodes).
	 */
	virtual inline void InitializeNodeInteraction( NodePointerType node );

	/** Initializing edges attributes. */
	virtual inline void InitializeEdge( EdgePointerType edge, const NodeType & source, const NodeType & target );

	/** Initializing nodes attributes. */
	virtual inline void InitializeNode( NodePointerType node );

	/** Updating a node's attribute. */
	virtual inline void UpdateNode( NodePointerType node );

	/** Initializing node weights without taking into account interations. */
	virtual inline void UpdateNodesWeightWithoutInteraction();

	/** Gets the fiber-like shape likelihood of a node (given its attributes). */
	virtual inline NodeWeightType GetFiberLikelihood( const NodeType & node ) const;

	/** Reset nodes connections. Correlates incoming and outgoing edges
	 * so they have the same index.
	 * It is called when the graph is being copied.
	 */
	virtual inline void ConnectNodeExtremities();

	/** Set outgoing edges attributes, from a node. */
	virtual inline void SetOutgoingEdgesWeight( const NodeType & node );

	/** Get the Weight of a node. */
	virtual inline NodeWeightType Weight( const NodeType & node ) const;

	/** Get the f attribute of an edge (between two nodes). */
	virtual inline double f( const NodeType & i, const NodeType & j ) const;

	/** Get the Force of a node. */
	virtual inline double Force( const NodeType & node ) const;

	/** Get the Force between two nodes. */
	virtual inline double Force( const NodeType & i, const NodeType & j, const EdgeType & ij ) const;

	/** Get geometric likelihood to be a fiber. */
	virtual inline double GeometricLikelihood( const NodeType & node ) const;

	/** Get the connection likelihood between two nodes. */
	virtual inline double ConnectionLikelihood( const NodeType & i, const NodeType & j, const EdgeType & edge ) const;

	/** Get the connection likelihood between two not connection nodes.
	 * Not implemented, returns 1.
	 */
	virtual inline double VirtualConnectionLikelihood( const NodeType & i, const NodeType & j ) const;

	/** Get the mean point P of a set of points between iterator begin and end. */
	template< class TIterator >
	inline VectorType P( TIterator begin, TIterator end ) const;

	/** Get the director vector V of a set of points between iterator begin and end. */
	template< class TIterator >
	inline VectorType V( TIterator begin, TIterator end, const VectorType & p ) const;

	/** Get the average distance of a set of point to a vector v passing through p. */
	template< class TIterator >
	inline double AverageDistance( TIterator begin, TIterator end, const VectorType & p, const VectorType & v ) const;

	/** Get the mean grayscale value of a set of pixel between begin and end. */
	template< class TIterator >
	inline double Mean( TIterator begin, TIterator end ) const;

	/** Get the average distance of the bresenham line from the director vector. */
	template< class TIterator >
	inline double AverageDistanceVariance( const VectorType & point, const VectorType & vector, double AverageDistance, TIterator begin, TIterator end ) const;


	/** Fill gaps if two lines are separated. Not tested yet, and not demonstrated
	 * usefull yet.
	 */
	inline void FillGaps( double threshold );

	/** Get the maximum node (considering its weight) below a threshold. */
	inline NodePointerType GetMaximumNodeWeight( const NodeWeightType & threshold );

	/** Merge indexes of nodes i an j. Indexes of i and j are copied */
	inline IndexContainerType MergeIndexes( const NodeType & i, const NodeType & j, bool ifront, bool jfront ) const;

	/** Get the minimum distance of index with the line (p,v). */
	inline double Distance( const IndexType & index, const VectorType & p, const VectorType & v ) const;

	/** Get the average distance of the set of points generated by the vector v,
	 * using bresenham algorithm.
	 */
	inline double TheoricalAverageDistance( const VectorType & v, double norm ) const;

	/** Tells if i is connected to j by the front side. */
	inline bool IsConnectedToFront( const NodeType & i, const NodeType & j ) const;

	/** Tells if i is connected to j by the back side. */
	inline bool IsConnectedToBack( const NodeType & i, const NodeType & j ) const;

	/** Tells if a node i is likely to be connected by the front side with j, even
	 * if i and j are not connected.
	 */
	inline bool IsLikelyConnectedToFront( const NodeType & i, const NodeType & j ) const;

	/** Tells if a node i is likely to be connected by the back side with j, even
	 * if i and j are not connected.
	 */
	inline bool IsLikelyConnectedToBack( const NodeType & i, const NodeType & j ) const;

	/** Tells if i is only connected to j by the front side. It can be connected
	 * to both sides at the same time is j is a loop.
	 */
	inline bool IsOnlyConnectedToFront( const NodeType & i, const NodeType & j ) const;

	/** Tells if i is only connected to j by the back side. It can be connected
	 * to both sides at the same time is j is a loop.
	 */
	inline bool IsOnlyConnectedToBack( const NodeType & i, const NodeType & j ) const;

	/** Tells if i is connected to itself (a loop). */
	inline bool IsSelfConnected( const NodeType & i ) const;

	/** Tells if i and j have a connection. */
	inline bool AreConnected( const NodeType & i, const NodeType & j ) const;

	/** Synchronize incoming an outgoing edges array of a node, so they have the same index. */
	inline void SynchronizeIncomingOnOutgoingEdges( NodePointerType node );

	/** Deletes every edges in both sides of node i, front or back side. */
	inline void DeleteLinks( NodePointerType & i, bool front );

	/** Returns the maximum distance of a set of point to a line defined by a point and a vector.
	 * the position of the point is set in pos. */
	inline double MaxDistance( const IndexContainerType & indexes, const VectorType & point, const VectorType & vector, unsigned int & pos ) const;

	/** Tells if two nodes should be merges given a spatial threshold. */
	inline bool CanBeMerged( const NodeType & source, const NodeType & target, double threshold ) const;

	/** Tells if two nodes should be split, given a spatial threshold. */
	inline bool CanBeSplit( const IndexContainerType & indexes, const VectorType & point, const VectorType & vector, double threshold, unsigned int & pos ) const;

	/* Connect nodes i and j. Has to be debugged! */
	inline void Connect( NodePointerType & i, NodePointerType & j, bool iFront, bool jFront );

	/** Tries to merge a node with one of its neighbors, given a threshold. */
	inline bool MergeSuccess( NodePointerType node, const double & threshold );

	/** Merge two nodes i and j. */
	inline bool Merge( NodePointerType & i, NodePointerType & j, bool ijOutgoingSide, bool jiOutgoingSide, bool reloopingPossible = true );

	/** Returns the maximum intersection label (FrontLabel and BackLabel of the nodes. */
	inline NodeIdentifierType MaxIntersectionLabel() const;

	/** Tries to split a line into two lines, if it is a broken line (given a spatial threshold. */
	NodeIdentifierType TryToSplit( NodePointerType node, double threshold, NodeIdentifierType & intersectionLabel );

	/** Splits a node into two. */
	inline NodePointerType Split( NodePointerType node, unsigned int pos, NodeIdentifierType inter );



protected:

	SkeletonGraph(){};
	~SkeletonGraph(){};


private:
	SkeletonGraph( const Self & ); // purposely not implemented
	void operator=( const Self & ); // purposely not implemeted

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSkeletonGraph.hxx"
#endif

#endif
