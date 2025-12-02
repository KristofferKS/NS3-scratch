set(target_prefix scratch_)

function(create_scratch source_files)
  # Return early if no sources in the subdirectory
  list(LENGTH source_files number_sources)
  if(number_sources EQUAL 0)
    return()
  endif()

  # If the scratch has more than a source file, we need to find the source with
  # the main function
  set(scratch_src)
  foreach(source_file ${source_files})
    file(READ ${source_file} source_file_contents)
    string(REGEX MATCHALL "main[(| (]" main_position "${source_file_contents}")
    if(CMAKE_MATCH_0)
      list(APPEND scratch_src ${source_file})
    endif()
  endforeach()

  list(LENGTH scratch_src scratch_src_len)

  # If there is no main function, raise an error
  if(${scratch_src_len} EQUAL 0)
    message(FATAL_ERROR "The following scratch source files do not contain a main function: ${source_files}")
  endif()

  # If there are multiple main functions, raise an error
  if(${scratch_src_len} GREATER 1)
    message(FATAL_ERROR "The following scratch source files contain ${scratch_src_len} files with a main function: ${scratch_src}")
  endif()

  # If there is a single main function, continue normally

  # Get parent directory name
  get_filename_component(scratch_dirname ${scratch_src} DIRECTORY)
  string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "" scratch_dirname
                 "${scratch_dirname}"
  )
  string(REPLACE "/" "_" scratch_dirname "${scratch_dirname}")

  # Get source name
  get_filename_component(scratch_name ${scratch_src} NAME_WLE)

  set(target_prefix scratch_)
  if(scratch_dirname)
    # Join the names together if dirname is not the scratch folder
    set(target_prefix scratch${scratch_dirname}_)
  endif()

  # Get source absolute path and transform into relative path
  get_filename_component(scratch_src ${scratch_src} ABSOLUTE)
  get_filename_component(scratch_absolute_directory ${scratch_src} DIRECTORY)
  string(REPLACE "${PROJECT_SOURCE_DIR}" "${CMAKE_OUTPUT_DIRECTORY}"
                 scratch_directory ${scratch_absolute_directory}
  )
  build_exec(
          EXECNAME ${scratch_name}
          EXECNAME_PREFIX ${target_prefix}
          SOURCE_FILES "${source_files}"
          LIBRARIES_TO_LINK "${ns3-libs}" "${ns3-contrib-libs}"
          EXECUTABLE_DIRECTORY_PATH ${scratch_directory}/
  )
endfunction()

# Scan *.cc files in ns-3-dev/scratch and build a target for each
file(GLOB single_source_file_scratches CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/[^.]*.cc)
foreach(scratch_src ${single_source_file_scratches})
  create_scratch(${scratch_src})
endforeach()

# --- Optional Python scratch helpers -------------------------------------------------
# Add custom targets for Python scratch scripts so they appear in the CMake Tools
# outline. They are convenience run targets; debugging still uses the VS Code
# "Python: Debug current ns-3 script" configuration. Adjust or add more scripts as needed.

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
  # tester.py run target
  add_custom_target(py_tester
    COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${CMAKE_OUTPUT_DIRECTORY}/bindings/python:$ENV{PYTHONPATH}
            LD_LIBRARY_PATH=${CMAKE_OUTPUT_DIRECTORY}/lib:$ENV{LD_LIBRARY_PATH}
            ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tester.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Run ns-3 Python scratch tester.py"
    VERBATIM
  )
endif()

# Auto-add custom run targets for every top-level Python script (*.py) in this
# scratch directory (excluding files starting with a dot). This lets CMake Tools
# show them in the target outline. They are convenience run targets only.
file(GLOB scratch_python_scripts CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/[^.]*.py)
foreach(pyfile ${scratch_python_scripts})
  get_filename_component(pyname ${pyfile} NAME_WE)
  string(REPLACE "." "_" pyname ${pyname})
  # Skip if we already defined an explicit target with same name
  if(TARGET py_${pyname})
    continue()
  endif()
  add_custom_target(py_${pyname}
    COMMAND ${CMAKE_COMMAND} -E env
            PYTHONPATH=${CMAKE_OUTPUT_DIRECTORY}/bindings/python:$ENV{PYTHONPATH}
            LD_LIBRARY_PATH=${CMAKE_OUTPUT_DIRECTORY}/lib:$ENV{LD_LIBRARY_PATH}
            ${Python3_EXECUTABLE} ${pyfile}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Run ns-3 Python scratch ${pyfile}"
    VERBATIM
  )
endforeach()

# Scan ns-3-dev/scratch subdirectories
file(
  GLOB scratch_subdirectories
  CONFIGURE_DEPENDS
  LIST_DIRECTORIES true
  ${CMAKE_CURRENT_SOURCE_DIR}/**
)
# Filter out files
foreach(entry ${scratch_subdirectories})
  if(NOT (IS_DIRECTORY ${entry}))
    list(REMOVE_ITEM scratch_subdirectories ${entry})
  endif()
endforeach()

# Build scratches per directory or following CMakeLists.txt instructions
foreach(subdir ${scratch_subdirectories})
  if(EXISTS ${subdir}/CMakeLists.txt)
    # If the subdirectory contains a CMakeLists.txt file
    # we let the CMake file manage the source files
    #
    # Use this if you want to link to external libraries
    # without creating a module
    add_subdirectory(${subdir})
  else()
    # Otherwise we pick all the files in the subdirectory
    # and create a scratch for them automatically
    file(GLOB scratch_sources CONFIGURE_DEPENDS ${subdir}/[^.]*.cc)
    create_scratch("${scratch_sources}")
  endif()
endforeach()
