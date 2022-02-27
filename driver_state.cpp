#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    
    state.image_color = new pixel[width * height];
    for (int i = 0; i < width*height; i++)
	state.image_color[i] = make_pixel(0,0,0);

    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{

    data_geometry* t = new data_geometry[3];
    auto ptr = state.vertex_data;

    data_vertex in{};

    switch (type) {

	case render_type::triangle:
		for ( int i = 0, j = 0; i < state.num_vertices; i++, j++ ) {
			t[j].data = ptr;
			in.data = ptr;
			state.vertex_shader(in, t[j], state.uniform_data);
			if ( j == 2 ) {
				rasterize_triangle(state, t[0], t[1], t[2]);
				j = -1;
			}
			ptr += state.floats_per_vertex;
		}
		break;
	case render_type::indexed:
		break;
	case render_type::fan:
		break;
	case render_type::strip:
		break;
	default:
		break;
    }
   
    delete [] t;
    //std::cout<<"TODO: implement rendering."<<std::endl;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state, v0, v1, v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{

    data_geometry* v = new data_geometry[3];
    int x[3], y[3];

    v[0] = v0;
    v[1] = v1;
    v[2] = v2;

    for ( int d = 0; d < 3; d++ ) {
	int i = (state.image_width / 2.0) * v[d].gl_Position[0] + (state.image_width / 2.0);
	int j = (state.image_height / 2.0) * v[d].gl_Position[1] + (state.image_height / 2.0);
	x[d] = i;
	y[d] = j;
	state.image_color[i+j*state.image_width] = make_pixel(255,255,255);
    }

    auto *data = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment fragData{data};

    auto output = new data_output;
/*
    for(int i = 0; i < state.floats_per_vertex; i++) {
	switch(state.interp_rules[i]) {
		case interp_type::flat:
			fragData.data[i] = v0.data[i];
			break;
		case interp_type::smooth:
			break;
		case interp_type::noperspective:
			break;
		default:
			break;
	}
    }

    state.fragment_shader(fragData, *output, state.uniform_data);
*/
    float areaABC = (0.5f * ((x[1]*y[2] - x[2]*y[1]) + (x[2]*y[0] - x[0]*y[2]) + (x[0]*y[1] - x[1]*y[0])));

    for(int j = 0; j < state.image_height; j++) {
        for(int i = 0; i < state.image_width; i++) {
            float alpha = (0.5f * ((x[1] * y[2] - x[2] * y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)) / areaABC;
            float beta =  (0.5f * ((x[2] * y[0] - x[0] * y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j)) / areaABC;
            float gamma = (0.5f * ((x[0] * y[1] - x[1] * y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j)) / areaABC;

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
		for(int k = 0; k < state.floats_per_vertex; k++) {
			switch(state.interp_rules[k]) {
				case interp_type::flat:
					fragData.data[k] = v0.data[k];
					break;
				case interp_type::smooth:
					break;
				case interp_type::noperspective:
					fragData.data[k] = (alpha * v0.data[k]) + (beta * v1.data[k]) + (gamma * v2.data[k]);
					break;
				default:
					break;
			}
	        }
	     
	     	state.fragment_shader(fragData, *output, state.uniform_data);
             	state.image_color[i + j * state.image_width] = make_pixel(output->output_color[0] * 255, output->output_color[1] * 255, output->output_color[2] * 255);

             }
    	}
    }

    delete [] data;
    delete output;
    delete [] v;


    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

