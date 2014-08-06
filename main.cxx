// THIS PROJECT USES C++11
// NO LOWER VERSION
//


#include <iostream>
#include <fstream>
#include <math.h>
#include <list>
#include <stack>


	/*
	 * Image
	 */

#include "itkImage.h"


const unsigned int Dim = 3;
typedef float RealType;

typedef unsigned short 					PixelType;
typedef itk::Image< PixelType, Dim > 	ImageType;
typedef itk::Image< RealType, Dim > 	RealImageType;
typedef itk::ImageRegion< Dim >			ImageRegionType;

const PixelType MaxPixelValue = itk::NumericTraits< PixelType >::max();
const PixelType MinPixelValue = itk::NumericTraits< PixelType >::min();

#include "itkSkeletonGraph.h"
#include "itkLineGraphTraits.h"

typedef itk::LineGraphTraits< Dim >				GraphTraitsType;
typedef itk::SkeletonGraph< GraphTraitsType >	GraphType;

	/*
	 * Inputs and Outputs
	 */

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileReader< RealImageType > RealReaderType;
typedef itk::ImageFileWriter< ImageType > WriterType;

	/*
	 * Filters
	 */

#include "itkBinaryThinningImageFilter3D.h"
#include "itkSkeletonToGraphFilter.h"
#include "itkGraphToSkeletonImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFiberExtractionGraphFilter.h"
#include "itkLineShapeImageFilter.h"

#include "itkAVLIntervalTreeNode.h"
#include "itkAVLVectorIntervalTreeNode.h"
#include "itkHalfConvexHullTreeNode.h"
#include "itkHalfConvexHullTreeContainer.h"
#include "itkConvexHull.h"
#include "itkBlurred3DLineSegmentation.hxx"

typedef itk::Vector< int, 2 >											VectorType;
typedef itk::AVLIntervalTreeNode< VectorType >									AVLTreeNodeType;
typedef itk::HalfConvexHullTreeNode< VectorType >					HalfConvexHullTreeNodeType;
typedef itk::HalfConvexHullTreeContainer< HalfConvexHullTreeNodeType >			HalfConvexHullTreeContainerType;
typedef itk::BinaryTreeContainer< AVLTreeNodeType >						BinaryTreeContainerType;
typedef AVLTreeNodeType													TreeNodeType;

typedef itk::LineShapeImageFilter< ImageType > 							LineShapeFilterType;
typedef itk::BinaryThinningImageFilter3D< ImageType, ImageType >		BinaryThinningFilterType;

typedef itk::ImageFileWriter< RealImageType >							RealWriterType;

typedef itk::SkeletonToGraphFilter< ImageType, GraphType, RealImageType, ImageType >	SkeletonToGraphFilterType;
typedef itk::GraphToSkeletonImageFilter< GraphType, RealImageType >			GraphToSkeletonFilterType;
typedef itk::FiberExtractionGraphFilter< GraphType >					FiberExtractionFilterType;

typedef itk::MultiplyImageFilter< ImageType, RealImageType, RealImageType > MultiplyFilterType;

void histoGraphMap( GraphType::Pointer graph );

void display( TreeNodeType* node );


int main( int argc, char* argv[] )
{
	/*
	 * Checking input parameters
	 * The parameter is the file without the extension
	 */

	if( argc < 2 )
	{
		std::cout << "Error : Not enough input arguments" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str = std::string( argv[1] );

	WriterType::Pointer writer = WriterType::New();
	RealWriterType::Pointer realWriter = RealWriterType::New();

	/*
	 *	Reading input image
	 */

	// The file format is nrrd, but it can be changed
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( str + ".nrrd" );
	reader->Update();


	std::cout << " ***** LineShapeFilter " << std::endl;
	// This filter does a hessian + gaussian filter
	// We have then access to eigen values

	LineShapeFilterType::Pointer lineShapeFilter = LineShapeFilterType::New();

	lineShapeFilter->SetInput( reader->GetOutput() );
	// this parameter has to be about the width of a fiber
	lineShapeFilter->SetSigma( 2.0f );
	lineShapeFilter->SetExtractBrightLine( false );
	lineShapeFilter->EigenValuesExtractionOn();
	lineShapeFilter->LabelImageOn();
	lineShapeFilter->Update();

	// Binary image of the fiber-like shapes
	writer->SetInput( lineShapeFilter->GetBinaryOutput() );
	writer->SetFileName( str + "_binary.nrrd" );
	writer->Update();

	// Label image depending of the sign of the 3 sorted eigenvalues
	writer->SetInput( lineShapeFilter->GetLabelOutput() );
	writer->SetFileName( str + "_label.nrrd" );
	writer->Update();

	// Eigen value 0 in grayscale maske with the binary output
	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 0 ) );
	realWriter->SetFileName( str + "_eigen0.nrrd" );
	realWriter->Update();

	// Eigen value 1 in grayscale maske with the binary output
	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 1 ) );
	realWriter->SetFileName( str + "_eigen1.nrrd" );
	realWriter->Update();

	// Eigen value 2 in grayscale maske with the binary output
	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 2 ) );
	realWriter->SetFileName( str + "_eigen2.nrrd" );
	realWriter->Update();
	
	// Vesselness output (fiber-shape likelyhood)
	realWriter->SetInput( lineShapeFilter->GetVesselnessOutput() );
	realWriter->SetFileName( str + "_vesselness.nrrd" );
	realWriter->Update();

	std::cout << " ***** Skeletonization " << std::endl;
	// We give the binary output of the last filter

	BinaryThinningFilterType::Pointer thinningFilter = BinaryThinningFilterType::New();
	thinningFilter->SetInput( lineShapeFilter->GetOutput() );
	thinningFilter->Update();

	writer->SetInput( thinningFilter->GetOutput() );
	writer->SetFileName( str + "_skeleton.nrrd" );
	writer->Update();

	std::cout << " ***** Vesselness Skeleton " << std::endl;
	// We get an image of the skeleton masked with the vesselness image

	MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
	multiplyFilter->SetInput1( thinningFilter->GetOutput() );
	multiplyFilter->SetInput2( lineShapeFilter->GetVesselnessOutput() );
	multiplyFilter->Update();

	realWriter->SetInput( multiplyFilter->GetOutput() );
	realWriter->SetFileName( str + "_skeleton_vesselness.nrrd" );
	realWriter->Update();

	std::cout << " ***** Skeleton To Graph Filter " << std::endl;
	// From the skeleton, we build a graph
	
	SkeletonToGraphFilterType::Pointer skeletonToGraphFilter = SkeletonToGraphFilterType::New();
	skeletonToGraphFilter->SetInput( thinningFilter->GetOutput() );
	skeletonToGraphFilter->SetVesselnessInput( lineShapeFilter->GetVesselnessOutput() );
	skeletonToGraphFilter->Update();

	// Identifier output of the skeleton to graph
	// This is a skeleton on which values are the number of the nodes
	itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::Pointer iwriter = itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::New();
	iwriter->SetInput( skeletonToGraphFilter->GetIdentifierOutput() );
	iwriter->SetFileName( str + "_identifier.nrrd" );
	iwriter->Update();

	std::cout << " ***** Fiber extraction Graph Filter  " << std::endl;
	// It takes the graph as input and processes it
	// This version works with straight fibers

	FiberExtractionFilterType::Pointer fiberExtractionFilter = FiberExtractionFilterType::New();
	fiberExtractionFilter->SetInput( skeletonToGraphFilter->GetOutput() );
	fiberExtractionFilter->Update();

	std::cout << " ***** Converting graph to image " << std::endl;
	// To see what we get
	// If you want to see any values, go into the GraphToSkeletonFilter
	// file, and change what you write (very small class).

	GraphToSkeletonFilterType::Pointer graphToSkeletonFilter = GraphToSkeletonFilterType::New();
	graphToSkeletonFilter->SetInput( fiberExtractionFilter->GetOutput() );
	graphToSkeletonFilter->SetRegion( reader->GetOutput()->GetLargestPossibleRegion() );
	graphToSkeletonFilter->Update();

	// Writing what we want to see
	realWriter->SetInput( graphToSkeletonFilter->GetOutput() );
	realWriter->SetFileName( str + "_fiber_extraction.nrrd" );
	realWriter->Update();

	return EXIT_SUCCESS;
}
