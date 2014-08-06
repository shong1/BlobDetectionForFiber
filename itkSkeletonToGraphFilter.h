#ifndef __itkSkeletonToGraphFilter_h
#define __itkSkeletonToGraphFilter_h

#include "Graph/itkGraphSource.h"
#include "itkObjectFactory.h"

#include "itkFixedArray.h"

namespace itk
{
/** \class SkeletonToGraphFilter
 * \brief Computes a graph from a binary skeleton input
 *
 * The filter takes as an input a binary skeleton image and a grayscale likelyhood image.
 * Each node of the graph represents a branch of the skeleton, connected with
 * other branches by edges. An index container stores all the indexes of the line.
 * Beside that, a node stores any wanted information calculated from the line and its connections.
 *
 * The input image type can be any binary image. The vesselness image type has to be
 * a grayscale image type, and the label image type should have integer type.
 *
 * The label image is an output image the labels every branch of the skeleton
 * with its corresponding node in the graph.
 *
 * The graph type has to be inhereted from GraphSource, a class from the Graph MIDAS library :
 * https://github.com/midas-journal/midas-journal-306/tree/master/Source
 *
 */
template< class TInputImage, class TOutputGraph, class TVesselnessImage, class TLabelImage = TInputImage >
class ITK_EXPORT SkeletonToGraphFilter : public GraphSource< TOutputGraph >
{
public:
	/** Standard class typedefs. */
	typedef SkeletonToGraphFilter			Self;
	typedef GraphSource< TOutputGraph >		Superclass;
	typedef SmartPointer< Self >			Pointer;
	typedef SmartPointer< const Self >		ConstPointer;

	/** Method for creation through the object factory. */
	itkNewMacro( Self );

	/** Run-time type information (and related methods). */
	itkTypeMacro( SkeletonToGraphFilter, GraphSource );

	/**  Some Image related typedefs. */
	typedef TInputImage							ImageType;
	typedef typename ImageType::Pointer			ImagePointer;
	typedef typename ImageType::ConstPointer	ImageConstPointer;
	typedef typename ImageType::RegionType		RegionType;
	typedef typename ImageType::PixelType		PixelType;
	typedef typename ImageType::IndexType		IndexType;
	typedef typename ImageType::OffsetType		OffsetType;

	/** Some Graph related typedefs. */
	typedef TOutputGraph								GraphType;
	typedef typename GraphType::GraphTraitsType			GraphTraitsType;
	typedef typename GraphType::Pointer					GraphPointer;
	typedef typename GraphType::NodeType				NodeType;
	typedef typename GraphType::NodePointerType			NodePointerType;
	typedef typename GraphType::NodeIdentifierType		NodeIdentifierType;
	typedef typename GraphTraitsType::NodeWeightType	NodeWeightType;
	typedef typename GraphType::EdgeType				EdgeType;
	typedef typename GraphType::EdgePointerType			EdgePointerType;
	typedef typename GraphType::EdgeIdentifierType		EdgeIdentifierType;
	typedef typename GraphTraitsType::EdgeWeightType	EdgeWeightType;
	typedef typename GraphTraitsType::VectorType		VectorType;
	typedef typename GraphTraitsType::PointType			PointType;
	typedef typename GraphTraitsType::IndexContainerType	IndexContainerType;
	typedef EdgeIdentifierType							IdentifierType;

	/** Some Label Image related typedefs. */
	typedef TLabelImage									LabelImageType;
	typedef typename LabelImageType::Pointer			LabelImagePointer;
	typedef typename LabelImageType::ConstPointer		LabelImageConstPointer;
	typedef typename LabelImageType::PixelType			LabelType;

	static const unsigned int ImageDimension = ImageType::ImageDimension;

	/** Some Identifier related typedefs. */
	typedef Image< NodeIdentifierType, ImageDimension >	IdentifierImageType;
	typedef typename IdentifierImageType::Pointer		IdentifierImagePointer;
	typedef typename IdentifierImageType::ConstPointer	IdentifierImageConstPointer;

	/** Some Vesselness related typedefs. */
	typedef TVesselnessImage							VesselnessImageType;
	typedef typename VesselnessImageType::Pointer		VesselnessImagePointer;
	typedef typename VesselnessImageType::ConstPointer	VesselnessImageConstPointer;

	/** Iterator typedefs */
	typedef ConstNeighborhoodIterator< ImageType >				ImageIteratorType;
	typedef NeighborhoodIterator< LabelImageType >				LabelImageIteratorType;
	typedef NeighborhoodIterator< IdentifierImageType >			IdentifierImageIteratorType;
	typedef ConstNeighborhoodIterator< VesselnessImageType >	VesselnessImageIteratorType;

	typedef std::vector< unsigned int > 				ArrayType;

	/** Set the input image of this process object. */
	void SetInput( unsigned int idx, const ImageType *input );
	void SetInput( ImageType * );

	/** Get the input image of this process object. */
	const ImageType * GetInput( unsigned int idx );
	
	/** Get the output graph of this process object. */
	GraphType * GetOutput( void );

	/** Get the output label image of this process object. */
	LabelImageType * GetLabelOutput( void );

	/** Get the output identifier image of this process object. */
	IdentifierImageType * GetIdentifierOutput( void );

	/** Set/Get the input vesselness image, which should be a grayscale image
	 * giving some likelihood.
	 */
	void SetVesselnessInput( const VesselnessImageType * vesselness );
	const VesselnessImageType * GetVesselnessInput() const;

	/** Prepare the output */
//	virtual void GenerateOutputInformation( void );

protected:
	SkeletonToGraphFilter();
	~SkeletonToGraphFilter();
//	void PrintSelf( std::ostream& os, Indent indent ) const;

	/** Generate data function. */
	virtual void GenerateData();

	/** Function to transform the skeleton image to a graph. ALL iterators shall
	 * be at the same location in the image, anywhere in the skeleton. It will
	 * transform the local branch to a node and connect it.
	 */
	void TransformSkeletonToGraph( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit );

	/** Function to make a node. It should be called when all the iterators are pointing
	 * the same pixel on the image. The node parameter should already be allocated
	 * using the New() routine. The boolean borderFound is here because the branch
	 * has to be went through once or twice. borderFound needs to be set to false at the
	 * first call of the function If the pixel is on the middle of the branch,
	 * MakeNode will go through one direction, and then the other. In some case, the
	 * function doesn't do the come back, and then borderFound is still set to false.
	 * if borderFound is still equals to false at the end of the function, MakeNode
	 * needs to be recalled, with borderFound = true, at the same location for the
	 * iterators as at the last call.
	 */
	void MakeNode( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, bool * borderFound );

	/** Function to call when a pixel at an intersection is found (>2 connected
	 * neighbors). It fills the iit iterator at the location of every connected
	 * pixels to the intersection. It is needed to know if an intersection is
	 * virgin or not, and to then connect nodes.
	 */
	void FillIntersection( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, LabelType & label );

	/** Function that will connect nodes at a non virgin intersection. If the
	 * intersection is non virgin, it will connect every connected nodes that
	 * are already created. It shall be called several times for one intersection.
	 * The label parameter is the label of the intersection. 
	 */
	void ConnectNodes( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit, NodePointerType node, LabelType label );

	/** Function to connect the node with itself. */
	void ConnectSelf( NodePointerType node );

	/** Get number of connected neighbor pixels to the pixel at position it. */
	unsigned int GetNeighborhoodCount( const ImageIteratorType & it ) const;

	/** Get a mapping of the location of all connected neighbor piexels to
	 * the location at pixel it. For example, if there are 2 neighbors, map
	 * would have a size of 2. map[0] is the index of the location of
	 * the first neighbor, map[1] of the segond.
	 */
	unsigned int GetNeighborhoodMap( const ImageIteratorType & it, ArrayType & map ) const;


	/** A fonction to call from the beginning iterator of the image.
	 * It goes through all the image, and calls TransformSkeletonToGraph
	 * everytime it meets a white pixel.
	 */
	bool ReachBranchSkeleton( ImageIteratorType & it, LabelImageIteratorType & lit, IdentifierImageIteratorType & iit, VesselnessImageIteratorType & vit ) const;

	DataObject::Pointer MakeOutput( unsigned int idx );




private:
	SkeletonToGraphFilter( const SkeletonToGraphFilter & ); // purposely not implemented
	void operator=( const SkeletonToGraphFilter & );		// purposely not implemented
	
	unsigned int m_NeighborhoodCount;
	LabelType m_Label;
	RegionType m_Region;
};

} // namesapce itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSkeletonToGraphFilter.hxx"
#endif

#endif
