# Dear ImGui
add_library(imgui STATIC 
	imgui.h
	imgui.cpp
	imgui_widgets.cpp
	imgui_demo.cpp
	imgui_draw.cpp
	imgui_tables.cpp
	backends/imgui_impl_glfw.h
	backends/imgui_impl_glfw.cpp
	backends/imgui_impl_opengl3.h
	backends/imgui_impl_opengl3.cpp
)

target_include_directories(imgui PUBLIC . backends)
#add_definitions(-DIMGUI_IMPL_OPENGL_LOADER_GLAD)
target_link_libraries(imgui LINK_PRIVATE glad glfw)
#set_property(TARGET imgui PROPERTY FOLDER "External")



