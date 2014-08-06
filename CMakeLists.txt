CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT( concrete )

SET(CMAKE_AUTOMOC ON)

#Find ITK
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
	INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

FIND_PACKAGE(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config.in
${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
ADD_CUSTOM_TARGET(doc
${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile WORKING_DIRECTORY
${CMAKE_CURRENT_BINARY_DIR} COMMENT "Generating API documentation with Doxygen" VERBATIM)
ENDIF(DOXYGEN_FOUND)

OPTION(BUILD_DOCUMENTATION "build documentation" ON)
OPTION(BUILD_DOXYGEN "build doxygen" ON)
OPTION(ITK_BUILD_ALL_MODULE "build all module" ON)

IF(CMAKE_COMPILER_IS_GNUCXX)
	ADD_DEFINITIONS(-fpermissive -std=gnu++0x)
ENDIF()

#Find Qt
#FIND_PACKAGE(Qt4 REQUIRED)
#INCLUDE(${QT_USE_FILE})

# Add Qt files
#ADD_DEFINITIONS(${QT_DEFINITIONS})

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} )

#SET( My_sources main.cxx )
#SET( My_hearders )
#SET( My_forms )
#SET( My_RCC )

#QT4_ADD_RESOURCES(My_RCC_srcs ${My_RCC})
#QT4_WRAP_UI( My_forms_headers ${My_forms} )
#QT4_WRAP_CPP( My_headers_moc ${My_headers} )

SET( OtherHeaders itkBinaryThinningImageFilter3D.h itkLineShapeImageFilter.h Graph/itkGraphSource.h itkSkeletonToGraphFilter.h Graph/itkDefaultGraphTraits.h itkLineGraphTraits.h Graph/itkGraph.h Graph/itkGraphDuplicator.h itkGraphToSkeletonImageFilter.h Graph/itkGraphToImageFilter.h Graph/itkInPlaceGraphFilter.h itkSkeletonGraph.h itkFiberExtractionGraphFilter.h itkBinaryTreeContainer.h itkBinaryTreeNodeBase.h itkSearchBinaryTreeNode.h itkAVLTreeNode.h itkAVLIntervalTreeNode.h
	itkAVLVectorIntervalTreeNode.h
	itkHalfConvexHullTreeNode.h itkHalfConvexHullTreeContainer.h itkConvexHull.h itkBlurred3DLineSegmentation.h ) # All files not using Qt
SET( OtherSources main.cxx itkBinaryThinningImageFilter3D.txx itkLineShapeImageFilter.hxx Graph/itkGraphSource.txx Graph/itkGraph.txx Graph/itkGraphDuplicator.txx itkGraphToSkeletonImageFilter.hxx Graph/itkGraphToImageFilter.txx Graph/itkInPlaceGraphFilter.txx itkSkeletonGraph.hxx itkFiberExtractionGraphFilter.hxx itkBinaryTreeContainer.hxx itkBinaryTreeNodeBase.hxx itkSearchBinaryTreeNode.hxx itkAVLTreeNode.hxx itkAVLIntervalTreeNode.hxx itkAVLVectorIntervalTreeNode.hxx itkHalfConvexHullTreeNode.hxx
	itkHalfConvexHullTreeContainer.hxx itkConvexHull.hxx itkBlurred3DLineSegmentation.hxx )


#ADD_EXECUTABLE(skull_stripping ${My_sources} ${My_headers_moc} ${My_forms_headers} ${My_RCC_srcs} ${OtherHeaders} ${OtherSources} )
ADD_EXECUTABLE(concrete ${OtherHeaders} ${OtherSources} )

TARGET_LINK_LIBRARIES(concrete ${ITK_LIBRARIES} )
