// THIS PROJECT USES C++11
// NO LOWER VERSION
//
#include <iostream>
#include <fstream>
#include <math.h>
#include <list>
#include <stack>
#include <string>


using namespace std;
using std::ofstream;
using std::ifstream;


	/*
	 * Image
	 */

#include "itkImage.h"

const unsigned int Dim = 3;
typedef float RealType;
typedef unsigned char BinType;

typedef unsigned short 					PixelType;
typedef itk::Image< PixelType, Dim > 	ImageType;
typedef itk::Image< RealType, Dim > 	RealImageType;
typedef itk::ImageRegion< Dim >			ImageRegionType;
typedef itk::Image< BinType, Dim > 		LabelType;

typedef itk::Image< unsigned long , Dim > 	SegLabelType;


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

typedef itk::ImageFileWriter< LabelType > LabelWriterType;

typedef itk::ImageFileWriter<  SegLabelType > SegLabelWriterType;

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
#include "itkImageDuplicator.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkScalarImageKmeansImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkStatisticsImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkGrayscaleDilateImageFilter.h"


#include "itkGradientMagnitudeImageFilter.h"
#include "itkWatershedImageFilter.h"

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

typedef itk::ImageDuplicator< RealImageType >	DuplicatorType;



typedef itk::MinimumMaximumImageCalculator< ImageType >					MinMaxFilterType;
typedef itk::StatisticsImageFilter< ImageType > 						StaticFilterType;
typedef itk::ScalarImageKmeansImageFilter< ImageType, LabelType > 		KMeansType;
typedef itk::StatisticsImageFilter< SegLabelType > 						LabelStaticFilterType;


typedef itk::ConnectedComponentImageFilter< LabelType, ImageType >		CCAFilterType;
typedef itk::RelabelComponentImageFilter< ImageType, LabelType > 		RelabelFilterType;
typedef itk::BinaryThresholdImageFilter< ImageType, LabelType > 		BinaryFilterType;
typedef itk::BinaryBallStructuringElement< LabelType::PixelType, 3 >	StructureElementType;
typedef itk::BinaryErodeImageFilter< LabelType, LabelType, StructureElementType > 			BinErodeFilterType;
typedef itk::GrayscaleDilateImageFilter< ImageType, ImageType, StructureElementType >		GrayDilateFilterType;
typedef itk::SubtractImageFilter< ImageType, ImageType > 				SubtractorType;

typedef itk::GradientMagnitudeImageFilter< ImageType, RealImageType > 	GradientMagnitudeFilterType;
typedef itk::WatershedImageFilter< RealImageType > 						WatershedFilterType;

//#include "vtkSmartPointer.h"
//#include "vtkImageData.h"
//#include "vtkPointData.h"
//#include "vtkFloatArray.h"
//#include "vtkPolyData.h"
//#include "vtkUnsignedShortArray.h"
//#include "vtkDataObject.h"
//#include "vtkContourFilter.h"
//#include "vtkLookupTable.h"
//#include "vtkImageMapToColors.h"
//#include "vtkProbeFilter.h"
//#include "vtkSmoothPolyDataFilter.h"
//#include "vtkPolyDataMapper.h"
//#include "vtkScalarBarActor.h"
//#include "vtkTextProperty.h"
//#include "vtkTextMapper.h"
//#include "vtkActor2D.h"
//#include "vtkActor.h"
//#include "vtkRenderer.h"
//#include "vtkRenderWindow.h"
//#include "vtkRenderWindowInteractor.h"
//#include "vtkStructuredPointsReader.h"
//
//#define VTK_CREATE( type, var ) vtkSmartPointer< type > var = vtkSmartPointer< type >::New()


void histoGraphMap( GraphType::Pointer graph );

void display( TreeNodeType* node );

int main( int argc, char* argv[] )
{
	/*
	 * Checking input parameters
	 * The parameter is the file without the extension
	 */

//	if( argc < 2 )
//	{
//		std::cout << "Error : Not enough input arguments" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	std::string str = std::string( argv[1] );
	std::string str = "/home/shong/Research/FibersInConcrete/Data/StiffFibers_cropped";

	WriterType::Pointer writer = WriterType::New();
	RealWriterType::Pointer realWriter = RealWriterType::New();
	LabelWriterType::Pointer labelWriter = LabelWriterType::New();
	SegLabelWriterType::Pointer segLabelWriter = SegLabelWriterType::New();

	/*
	 *	Reading input image
	 */

	// The file format is nrrd, but it can be changed
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( str + ".nrrd" );
	reader->Update();
	ImageType::Pointer input = reader->GetOutput();

	StaticFilterType::Pointer statFilter = StaticFilterType::New();
	statFilter->SetInput( input );
	statFilter->Update();

	float meanI = statFilter->GetMean();
	float stdI = statFilter->GetSigma();
	float minImg = statFilter->GetMinimum();
	float maxImg = statFilter->GetMaximum();
	cout << "Mean : " << meanI << " Std : " << stdI << " Min : " << minImg << " Max : " << maxImg << endl;

	float statMin = meanI - 3 * stdI;
	float statMax = meanI + 3 * stdI;
	cout << "StatMin : " << statMin << " StatMax : " << statMax << endl;

	ImageType::RegionType region = input->GetLargestPossibleRegion();
	ImageType::SizeType dim3 = region.GetSize();

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		ImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;


		unsigned short val = input->GetPixel( idx );

		if( val < statMin )
			input->SetPixel( idx, statMin );
		else if( val > statMax )
			input->SetPixel( idx, statMax );
	}

	input->Update();

	writer->SetInput( input );
	writer->SetFileName( str + "_IntensityCut.nrrd" );
	writer->Update();


	GradientMagnitudeFilterType::Pointer gradientFilter = GradientMagnitudeFilterType::New();
	gradientFilter->SetInput( input );
	gradientFilter->Update();

//	realWriter->SetInput( gradientFilter->GetOutput() );
//	realWriter->SetFileName( str + "_gradMag.nrrd" );
//	realWriter->Update();

	WatershedFilterType::Pointer watersheder = WatershedFilterType::New();
	watersheder->SetInput( gradientFilter->GetOutput() );
	watersheder->SetThreshold( 0.001 );
	watersheder->SetLevel( 0.2 );
	watersheder->Update();

	segLabelWriter->SetInput( watersheder->GetOutput() );
	segLabelWriter->SetFileName( str + "_watershed_001_2.nrrd" );
	segLabelWriter->Update();

	LabelStaticFilterType::Pointer labelStatFilter = LabelStaticFilterType::New();
	labelStatFilter->SetInput( watersheder->GetOutput() );
	labelStatFilter->Update();

	SegLabelType::Pointer waterShedLabel = watersheder->GetOutput();


	float meanL = labelStatFilter->GetMean();
	float stdL = labelStatFilter->GetSigma();
	float minL = labelStatFilter->GetMinimum();
	float maxL = labelStatFilter->GetMaximum();

	cout << "Mean : " << meanL << " Std : " << stdL << " Min : " << minL << " Max : " << maxL << endl;

	double* meanArr = new double [ (int)maxL ];
	int* sizeArr = new int[ (int)maxL ];

	for( int i = 1; i < maxL; i++ )
	{
		meanArr[ i ] = 0;
		sizeArr[ i ] = 0;
	}

	for( int x = 0; x < dim3[0]; x++ )
	for( int y = 0; y < dim3[1]; y++ )
	for( int z = 0; z < dim3[2]; z++ )
	{
		ImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;

		SegLabelType::IndexType idxSeg;
		idxSeg[0] = x;
		idxSeg[1] = y;
		idxSeg[2] = z;

		unsigned short val = input->GetPixel( idx );
		unsigned long lVal = waterShedLabel->GetPixel( idxSeg );

		int arrIdx = (int)( lVal - 1 );

		double cMean = meanArr[ arrIdx ];
		int cSize = sizeArr[ arrIdx ] ;

		meanArr[ arrIdx ] = cMean + val;
		sizeArr[ arrIdx ] = cSize + 1;

//		if( x % 10 == 0 && y % 10 == 0 && z % 30 == 0 )
//		{
//			cout << "Array Value : " << meanArr[ arrIdx ] << " Input : " << val << endl;
//			cout << "Array Idx : " << arrIdx << " x, y, z : " << x << ", " << y << ", "  << z << endl;
//			cout << "Array Size : " << cSize << endl;
//		}
	}

//	string checkTxt = "/home/shong/Research/FibersInConcrete/check.csv";
//	ofstream out( checkTxt );

	for( int i = 0; i < maxL; i++ )
	{
		double meanVal = meanArr[ i ];
		double sizeVal = (double) sizeArr[ i ] ;

		if
		meanArr[ i ] = meanVal / sizeVal;

//		if( i % 30 == 0 || i == 12087 )
//		{
//			double vval = meanArr[ i ];
//			int isize = sizeArr[ i ];
//
//			out << i << " : " << vval << " Size : " << isize << endl;
//
////			cout << "Array Idx : " << i << endl;
////			cout << "Array Value : " << vval << endl;
////			cout << "Size : " << isize  << endl;
//		}
	}
//
//	out.close();



//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.002 );
//	watersheder->SetLevel( 0.2 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_2.nrrd" );
//	segLabelWriter->Update();
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.003 );
//	watersheder->SetLevel( 0.3 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_3.nrrd" );
//	segLabelWriter->Update();
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.004 );
//	watersheder->SetLevel( 0.4 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_4.nrrd" );
//	segLabelWriter->Update();
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.005 );
//	watersheder->SetLevel( 0.5 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_5.nrrd" );
//	segLabelWriter->Update();
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.006 );
//	watersheder->SetLevel( 0.6 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_6.nrrd" );
//	segLabelWriter->Update();
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.007 );
//	watersheder->SetLevel( 0.7 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_7.nrrd" );
//	segLabelWriter->Update();
//
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.008 );
//	watersheder->SetLevel( 0.8 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_8.nrrd" );
//	segLabelWriter->Update();
//
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.009 );
//	watersheder->SetLevel( 0.9 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_9.nrrd" );
//	segLabelWriter->Update();
//
//	watersheder->SetInput( gradientFilter->GetOutput() );
//	watersheder->SetThreshold( 0.01 );
//	watersheder->SetLevel( 1.0 );
//	watersheder->Update();
//
//	segLabelWriter->SetInput( watersheder->GetOutput() );
//	segLabelWriter->SetFileName( str + "_watershed_10.nrrd" );
//	segLabelWriter->Update();

/*
	// Binarization -Thresholding
	BinaryFilterType::Pointer binaryFilter = BinaryFilterType::New();
	binaryFilter->SetInput( input );
	binaryFilter->SetInsideValue( 255 );
	binaryFilter->SetOutsideValue( 0 );

	BinaryFilterType::Pointer binaryFilterL = BinaryFilterType::New();
	binaryFilterL->SetInput( input );
	binaryFilterL->SetInsideValue( 255 );
	binaryFilterL->SetOutsideValue( 0 );

	BinaryFilterType::Pointer binaryFilterH = BinaryFilterType::New();
	binaryFilterH->SetInput( input );
	binaryFilterH->SetInsideValue( 255 );
	binaryFilterH->SetOutsideValue( 0 );

	// Connected Component Analysis for Extremal Region Extraction
	CCAFilterType::Pointer ccaFilter = CCAFilterType::New();
	RelabelFilterType::Pointer sizeFilter = RelabelFilterType::New();
	RelabelFilterType::ObjectSizeType  minSize = 5;
	sizeFilter->SetMinimumObjectSize( minSize );

	CCAFilterType::Pointer ccaFilterL = CCAFilterType::New();
	RelabelFilterType::Pointer sizeFilterL = RelabelFilterType::New();
	sizeFilterL->SetMinimumObjectSize( minSize );

	CCAFilterType::Pointer ccaFilterH = CCAFilterType::New();
	RelabelFilterType::Pointer sizeFilterH = RelabelFilterType::New();
	sizeFilterH->SetMinimumObjectSize( minSize );

	// Erosion Filter for Boundary Extraction
	StructureElementType erodeBall;
	erodeBall.SetRadius( 1 );
	erodeBall.CreateStructuringElement();

	BinErodeFilterType::Pointer erodeFilter = BinErodeFilterType::New();
	erodeFilter->SetKernel( erodeBall );
	BinErodeFilterType::Pointer erodeFilterL = BinErodeFilterType::New();
	erodeFilterL->SetKernel( erodeBall );
	BinErodeFilterType::Pointer erodeFilterH = BinErodeFilterType::New();
	erodeFilterH->SetKernel( erodeBall );


	// Erosion for CCA labels
	GrayDilateFilterType::Pointer lDilateFilter = GrayDilateFilterType::New();
	lDilateFilter->SetKernel( erodeBall );

	GrayDilateFilterType::Pointer lDilateFilterL = GrayDilateFilterType::New();
	lDilateFilterL->SetKernel( erodeBall );

	GrayDilateFilterType::Pointer lDilateFilterH = GrayDilateFilterType::New();
	lDilateFilterH->SetKernel( erodeBall );

	SubtractorType::Pointer subtractor = SubtractorType::New();

	int d = 100;

	for( int i = statMin; i < statMax; i += 50 )
	{
		stringstream ssIdx;
		ssIdx << i;
		string sIdx = ssIdx.str();

		int iL = i - d;
		int iH = i + d;

		// Find Qs
		// Binarization - Thresholding
		binaryFilter->SetLowerThreshold( 0 );
		binaryFilter->SetUpperThreshold( i );
		binaryFilter->Update();
		LabelType::Pointer binLab = binaryFilter->GetOutput();

		ccaFilter->SetInput( binLab );
		ccaFilter->Update();

		int nObj = ccaFilter->GetObjectCount();

		cout << "No. Object : " << nObj << endl;



//		if( nObj < 3 )
//			break;
//
//		ImageType::Pointer Q = ccaFilter->GetOutput();
//
//		writer->SetInput(Q);
//		writer->SetFileName( str + "_checkQ.nrrd" );
//		writer->Update();
//
//		// Get boundary
//		lDilateFilter->SetInput( Q );
//		lDilateFilter->Update();
//		ImageType::Pointer DilatedLab = lDilateFilter->GetOutput();
//
//		writer->SetInput(DilatedLab);
//		writer->SetFileName( str + "_checkQE.nrrd" );
//		writer->Update();
//
//		subtractor->SetInput1( DilatedLab );
//		subtractor->SetInput2( Q );
//		subtractor->Update();
//
//		ImageType::Pointer dQ = subtractor->GetOutput();
//
//		writer->SetInput(subtractor->GetOutput());
//		writer->SetFileName( str + "_checkdQ.nrrd" );
//		writer->Update();
//
//		// Find Extremal Regions at i
//
//
//		// Q_i+delta, Q_i-delta
//		binaryFilterL->SetLowerThreshold( 0 );
//		binaryFilterL->SetUpperThreshold( iL );
//		binaryFilterL->Update();
//		LabelType::Pointer binLabL = binaryFilter->GetOutput();
//
//		binaryFilterH->SetLowerThreshold( 0 );
//		binaryFilterH->SetUpperThreshold( iH );
//		binaryFilterH->Update();
//		LabelType::Pointer binLabH = binaryFilter->GetOutput();
//
//		ccaFilterL->SetInput( binLabL );
//		ccaFilterL->Update();
//		ccaFilterH->SetInput( binLabH );
//		ccaFilterH->Update();
//
//		int nObjL = ccaFilterL->GetObjectCount();
//		int nObjH = ccaFilterH->GetObjectCount();
//
//		ImageType::Pointer Q_L = ccaFilterL->GetOutput();
//		ImageType::Pointer Q_H = ccaFilterH->GetOutput();

//		lDilateFilterL->SetInput( Q_L );
//		lDilateFilterL->Update();
//		ImageType::Pointer DilatedLabdL = lDilateFilterL->GetOutput();
//
//		lDilateFilterH->SetInput( Q_H );
//		lDilateFilterH->Update();
//		ImageType::Pointer DilatedLabH = lDilateFilterH->GetOutput();

	}
*/

	cout << "Process Done" << endl;


//	std::cout << " ***** LineShapeFilter " << std::endl;
//
//	// This filter does a hessian + gaussian filter
//	// We have then access to eigen values
//	LineShapeFilterType::Pointer lineShapeFilter = LineShapeFilterType::New();
//
//	lineShapeFilter->SetInput( reader->GetOutput() );
//	// this parameter has to be about the width of a fiber
//	lineShapeFilter->SetSigma( 5.0f );
//	lineShapeFilter->SetExtractBrightLine( false );
//	lineShapeFilter->EigenValuesExtractionOn();
//	lineShapeFilter->LabelImageOn();
//	lineShapeFilter->Update();

//	// Binary image of the fiber-like shapes
//	writer->SetInput( lineShapeFilter->GetBinaryOutput() );
//	writer->SetFileName( str + "_binary.nrrd" );
//	writer->Update();
//
//	// Label image depending of the sign of the 3 sorted eigenvalues
//	writer->SetInput( lineShapeFilter->GetLabelOutput() );
//	writer->SetFileName( str + "_label.nrrd" );
//	writer->Update();
//
//	// Eigen value 0 in grayscale maske with the binary output
//	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 0 ) );
//	realWriter->SetFileName( str + "_eigen0.nrrd" );
//	realWriter->Update();
//
//	// Eigen value 1 in grayscale maske with the binary output
//	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 1 ) );
//	realWriter->SetFileName( str + "_eigen1.nrrd" );
//	realWriter->Update();
//
//	// Eigen value 2 in grayscale maske with the binary output
//	realWriter->SetInput( lineShapeFilter->GetEigenValuesOutput( 2 ) );
//	realWriter->SetFileName( str + "_eigen2.nrrd" );
//	realWriter->Update();



/*
	RealImageType::Pointer eigImg0 = lineShapeFilter->GetEigenValuesOutput( 0 );
	RealImageType::Pointer eigImg1 = lineShapeFilter->GetEigenValuesOutput( 1 );
	RealImageType::Pointer eigImg2 = lineShapeFilter->GetEigenValuesOutput( 2 );

	RealImageType::RegionType region = eigImg0->GetLargestPossibleRegion();
	RealImageType::SizeType dim = region.GetSize();

	RealImageType::Pointer blobImg = RealImageType::New();
	blobImg->SetRegions( eigImg0->GetLargestPossibleRegion() );
	blobImg->SetSpacing( eigImg0->GetSpacing() );
	blobImg->SetOrigin( eigImg0->GetOrigin() );
	blobImg->Allocate();

	LabelType::Pointer blobBin = LabelType::New();
	blobBin->SetRegions( eigImg0->GetLargestPossibleRegion() );
	blobBin->SetSpacing( eigImg0->GetSpacing() );
	blobBin->SetOrigin( eigImg0->GetOrigin() );
	blobBin->Allocate();

	RealImageType::Pointer vesselImg = RealImageType::New();
	vesselImg->SetRegions( eigImg0->GetLargestPossibleRegion() );
	vesselImg->SetSpacing( eigImg0->GetSpacing() );
	vesselImg->SetOrigin( eigImg0->GetOrigin() );
	vesselImg->Allocate();


	cout << dim[0] << ", " << dim[1] << ", " << dim[2] << endl;

	float alpha = 0.5;
	float beta = 0.5;
	float c = 0.5;

	for( int x = 0; x < dim[0]; x++ )
	for( int y = 0; y < dim[1]; y++ )
	for( int z = 0; z < dim[2]; z++ )
	{
		RealImageType::IndexType idx;
		idx[0] = x;
		idx[1] = y;
		idx[2] = z;
		float lambda1 = eigImg0->GetPixel( idx );
		float lambda2 = eigImg1->GetPixel( idx );
		float lambda3 = eigImg2->GetPixel( idx );

		// eigImg2 > eigImg1 > eigImg0
		// lambda3 > lambda2 > lambda1
		float blVal = abs( lambda1 ) / sqrt( abs( lambda2 * lambda3 ) );
		float sVal = sqrt( lambda1 * lambda1 + lambda2 * lambda2 + lambda3 * lambda3 );
		float aVal = abs( lambda2 / lambda3 );

		float V0 = exp( -( aVal*aVal ) / ( 2 * alpha * alpha ) ) * ( 1 - exp( -(blVal * blVal) / ( 2 * beta * beta ) ) ) * ( 1 - exp( -( sVal * sVal ) / ( 2 * c * c ) ) );

		float Ves = ( 1 - exp( -( aVal*aVal ) / ( 2 * alpha * alpha ) ) ) * ( exp( -(blVal * blVal) / ( 2 * beta * beta ) ) ) * ( 1 - exp( -( sVal * sVal ) / ( 2 * c * c ) ) );
//		float V0 = 1 - exp( ( -blVal * blVal ) / ( 2 * beta * beta ) );
		blobImg->SetPixel( idx, V0 );
		vesselImg->SetPixel( idx, blVal );

		if( V0 > 0.2 && lambda1 < 0 && lambda2 < 0 && lambda3 < 0 )
			blobBin->SetPixel( idx, 255 );
		else
			blobBin->SetPixel( idx, 0 );
	}

	blobImg->Update();
	blobBin->Update();
	vesselImg->Update();

	cout << "Blob Detected?" << endl;

	realWriter->SetInput( blobImg );
	realWriter->SetFileName( str + "blobness.nrrd" );
	realWriter->Update();

	realWriter->SetInput( vesselImg );
	realWriter->SetFileName( str + "vesselnessF.nrrd" );
	realWriter->Update();


	// Label image depending of the sign of the 3 sorted eigenvalues
//	writer->SetInput( blobBin );
//	writer->SetFileName( str + "_label.nrrd" );
//	writer->Update();

	cout << "Blob Written" << endl;

	CCAFilterType::Pointer ccaFilter = CCAFilterType::New();
	ccaFilter->SetInput( blobBin );
	ccaFilter->Update();

	cout << "# Objects : " << ccaFilter->GetObjectCount() << endl;

	RelabelFilterType::Pointer sizeFilter = RelabelFilterType::New();
	RelabelFilterType::ObjectSizeType  minSize = 10;
	sizeFilter->SetInput( ccaFilter->GetOutput () );
	sizeFilter->SetMinimumObjectSize( minSize );
	sizeFilter->Update();

	labelWriter->SetInput( sizeFilter->GetOutput() );
	labelWriter->SetFileName( str + "labelS.nrrd" );
	labelWriter->Update();

	// Vesselness output (fiber-shape likelihood)
	realWriter->SetInput( lineShapeFilter->GetVesselnessOutput() );
	realWriter->SetFileName( str + "_vesselness.nrrd" );
	realWriter->Update();

	realWriter->SetInput( blobImg );
	realWriter->SetFileName( str + "_blobness.nrrd" );
	realWriter->Update();


	// Label image depending of the sign of the 3 sorted eigenvalues
	writer->SetInput( lineShapeFilter->GetLabelOutput() );
	writer->SetFileName( str + "_label.nrrd" );
	writer->Update();

	cout << "Process Done"<< endl;
*/


//
//	KMeansType::Pointer KmeanFilter = KMeansType::New();
//	KmeanFilter->SetInput( input );
//	KmeanFilter->SetUseNonContiguousLabels( true );
//	KmeanFilter->AddClassWithInitialMean( statMin );
//	KmeanFilter->AddClassWithInitialMean( meanI );
//	KmeanFilter->AddClassWithInitialMean( statMax );
//	KmeanFilter->Update();
//
//	KMeansType::ParametersType estimatedMean = KmeanFilter->GetFinalMeans();
//
//	int nClass = estimatedMean.Size();
//
//	cout << "# of clusters : " << nClass << endl;
//
//	labelWriter->SetInput( KmeanFilter->GetOutput() );
//	labelWriter->SetFileName( str + "kmeansLabel.nrrd" );
//	labelWriter->Update();


//	/***********************************************
//	 * Blobs are a problem in the image. They split
//	 * fibers into pieces. Most of the blobs have
//	 * the same label as fibers, but are isolated.
//	 * If we can segment such blobs, we can dilate
//	 * them and only them so they touch the rest
//	 * of the skeleton they might have split.
//	 * If I is the binary image and B the binary
//	 * image with only the blobs, D the
//	 * dilatation operator, and K the skeletonization
//	 * operator, K( I or D(B) ) should get rid of
//	 * the split artifact.
//	 ***********************************
//	***********/
//
//
//	std::cout << " ***** Skeletonization " << std::endl;
//	// We give the binary output of the last filter
//
//	BinaryThinningFilterType::Pointer thinningFilter = BinaryThinningFilterType::New();
//	thinningFilter->SetInput( lineShapeFilter->GetOutput() );
//	thinningFilter->Update();
//
//	writer->SetInput( thinningFilter->GetOutput() );
//	writer->SetFileName( str + "_skeleton.nrrd" );
//	writer->Update();
//
//	std::cout << " ***** Vesselness Skeleton " << std::endl;
//	// We get an image of the skeleton masked with the vesselness image
//
//	MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
//	multiplyFilter->SetInput1( thinningFilter->GetOutput() );
//	multiplyFilter->SetInput2( lineShapeFilter->GetVesselnessOutput() );
//	multiplyFilter->Update();
//
//	realWriter->SetInput( multiplyFilter->GetOutput() );
//	realWriter->SetFileName( str + "_skeleton_vesselness.nrrd" );
//	realWriter->Update();
//
//	std::cout << " ***** Skeleton To Graph Filter " << std::endl;
//	// From the skeleton, we build a graph
//
//	SkeletonToGraphFilterType::Pointer skeletonToGraphFilter = SkeletonToGraphFilterType::New();
//	skeletonToGraphFilter->SetInput( thinningFilter->GetOutput() );
//	skeletonToGraphFilter->SetVesselnessInput( lineShapeFilter->GetVesselnessOutput() );
//	skeletonToGraphFilter->Update();
//
//	// Identifier output of the skeleton to graph
//	// This is a skeleton on which values are the number of the nodes
//	itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::Pointer iwriter = itk::ImageFileWriter< SkeletonToGraphFilterType::IdentifierImageType >::New();
//	iwriter->SetInput( skeletonToGraphFilter->GetIdentifierOutput() );
//	iwriter->SetFileName( str + "_identifier.nrrd" );
//	iwriter->Update();
//
//	std::cout << " ***** Fiber extraction Graph Filter  " << std::endl;
//	// It takes the graph as input and processes it
//	// This version works with straight fibers
//
//	FiberExtractionFilterType::Pointer fiberExtractionFilter = FiberExtractionFilterType::New();
//	fiberExtractionFilter->SetInput( skeletonToGraphFilter->GetOutput() );
//	fiberExtractionFilter->Update();
//
//	std::cout << " ***** Converting graph to image " << std::endl;
//	// To see what we get
//	// If you want to see any values, go into the GraphToSkeletonFilter
//	// file, and change what you write (very small class).
//
//	GraphToSkeletonFilterType::Pointer graphToSkeletonFilter = GraphToSkeletonFilterType::New();
//	graphToSkeletonFilter->SetInput( fiberExtractionFilter->GetOutput() );
//	graphToSkeletonFilter->SetRegion( reader->GetOutput()->GetLargestPossibleRegion() );
//	graphToSkeletonFilter->Update();
//
//	// Writing what we want to see
//	realWriter->SetInput( graphToSkeletonFilter->GetOutput() );
//	realWriter->SetFileName( str + "_fiber_extraction.nrrd" );
//	realWriter->Update();

	return EXIT_SUCCESS;
}
