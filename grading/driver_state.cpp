#include "driver_state.h"
#include <cstring>
#include <limits>

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
    state.image_depth = new float [width * height];
    for (int i = 0; i < width*height; i++) {
	state.image_color[i] = make_pixel(0,0,0);
	state.image_depth[i] = std::numeric_limits<float>::max();
    }
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
    data_geometry t[3];
    data_geometry tmpT[3];
    data_vertex in[3];

    switch (type) {
	case render_type::triangle:
	{
		int ptr = 0;
		for ( int i = 0; i < state.num_vertices / 3; i++) {
			for (int j = 0; j < 3; j++, ptr += state.floats_per_vertex) {
				in[j].data = &state.vertex_data[ptr];
				tmpT[j].data = in[j].data;
				state.vertex_shader(in[j], tmpT[j], state.uniform_data);
				t[j] = tmpT[j];
			}
			clip_triangle(state, t[0], t[1], t[2], 0);
		}
		break;
	}
	case render_type::indexed:
		for (int i = 0; i < 3 * state.num_triangles; i += 3) {
			for (int j = 0; j < 3; j++) {
				in[j].data = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
				tmpT[j].data = in[j].data;
				state.vertex_shader(in[j], tmpT[j], state.uniform_data);
				t[j] = tmpT[j];
			}
			clip_triangle(state, t[0], t[1], t[2], 0);
		}
		break;
	case render_type::fan:
		for (int i = 0; i < state.num_vertices; i++) {
			for (int j = 0; j < 3; j++) {
				if (j != 0) in[j].data = state.vertex_data + ((state.floats_per_vertex) * (i+j));
				else in[j].data = state.vertex_data + (j * state.floats_per_vertex);
			
				tmpT[j].data = in[j].data;
				state.vertex_shader(in[j], tmpT[j], state.uniform_data);
				t[j] = tmpT[j];
			}
			clip_triangle(state, t[0], t[1], t[2], 0);
		}
		break;
	case render_type::strip:
		for (int i = 0; i < state.num_vertices - 2; i++) {
			for (int j = 0; j < 3; j++) {
				in[j].data = &state.vertex_data[(i+j) * state.floats_per_vertex];
				tmpT[j].data = in[j].data;
				state.vertex_shader(in[j], tmpT[j], state.uniform_data);
				t[j] = tmpT[j];
			}
			clip_triangle(state, t[0], t[1], t[2], 0);
		}
		break;
	default:
		break;
    }
    
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
    else {
	bool isClip = false;
	float A1, B1, B2; vec4 P1, P2;

	data_geometry first[3]; data_geometry second[3];
	vec4 A = v0.gl_Position; vec4 B = v1.gl_Position; vec4 C = v2.gl_Position;

	if (A[2] < -A[3] && B[2] < -B[3] && C[2] < -C[3]) return;
	else {
		if (A[2] < -A[3] && B[2] >= -B[3] && C[2] >= -C[3]) {
			isClip = true;
			B1 = (-B[3] - B[2]) / (A[2] + A[3] - B[3] - B[2]);
			B2 = (-A[3] - A[2]) / (C[2] + C[3] - A[3] - A[2]);
			P1 = B1 * A + (1 - B1) * B;
			P2 = B2 * C + (1 - B2) * A;

			first[0].data = new float[state.floats_per_vertex];
			first[1] = v1;
			first[2] = v2;

			for(int i = 0; i < state.floats_per_vertex; i++) {
				switch (state.interp_rules[i]) {
					case interp_type::flat:
						first[0].data[i] = v0.data[i];
						break;
					case interp_type::smooth:
						first[0].data[i] = B2 * v2.data[i] + (1 - B2) * v0.data[i];
						break;
					case interp_type::noperspective:
						A1 = B2 * v2.gl_Position[3] / (B2 * v2.gl_Position[3] + (1 - B2) * v0.gl_Position[3]);
						first[0].data[i] = A1 * v2.data[i] + (1 - A1) * v0.data[i];
						break;
					default:
						break;
					}
			}
			first[0].gl_Position = P2;
			clip_triangle(state, first[0], first[1], first[2], face+1);


			second[0].data = new float[state.floats_per_vertex];
			second[1] = v1;
			second[2] = first[0];

                        for(int i = 0; i < state.floats_per_vertex; i++) {
                                switch (state.interp_rules[i]) {
                                        case interp_type::flat:
                                                second[0].data[i] = v0.data[i];
                                                break;
                                        case interp_type::smooth:
                                                second[0].data[i] = B1 * v0.data[i] + (1 - B2) * v1.data[i];
                                                break;
                                        case interp_type::noperspective:
                                                A1 = B1 * v0.gl_Position[3] / (B1 * v0.gl_Position[3] + (1 - B1) * v1.gl_Position[3]);
                                                second[0].data[i] = A1 * v0.data[i] + (1 - A1) * v1.data[i];
                                                break;
                                        default:
                                                break;
                                        }
			}
			second[0].gl_Position = P1;
			clip_triangle(state, second[0], second[1], second[2], face+1);
		}
		if (isClip == false ) clip_triangle(state, v0, v1, v2, face+1);
	}
    }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    data_geometry* v = new data_geometry[3];
    v[0] = v0; v[1] = v1; v[2] = v2;

    float x[3], y[3], z[3];

    for ( int d = 0; d < 3; d++ ) {
	float i = (state.image_width / 2) * (v[d].gl_Position[0] / v[d].gl_Position[3]) + ((state.image_width / 2) - 0.5f);
	float j = (state.image_height / 2) * (v[d].gl_Position[1] / v[d].gl_Position[3]) + ((state.image_height / 2) - 0.5f);
	float k = (state.image_width / 2) * (v[d].gl_Position[2] / v[d].gl_Position[3]) + ((state.image_width / 2) - 0.5f);
	x[d] = i;
	y[d] = j;
	z[d] = k;
    }

    float minX = std::min(std::min(x[0], x[1]), x[2]);
    float maxX = std::max(std::max(x[0], x[1]), x[2]);
    float minY = std::min(std::min(y[0], y[1]), y[2]);
    float maxY = std::max(std::max(y[0], y[1]), y[2]);

    if (minX < 0) minX = 0;
    if (minY < 0) minY = 0;
    if (maxX > state.image_width) maxX = state.image_width;
    if (maxY > state.image_height) maxY = state.image_height;


    float area = ( (x[1]*y[2]) - (x[2]*y[1]) ) + ( (x[2]*y[0]) - (x[0]*y[2]) )+ ( (x[0]*y[1]) - (x[1]*y[0]) );

    data_output output;
    data_fragment fragData;
    float data[MAX_FLOATS_PER_VERTEX];
    

    for(int j = minY; j < maxY; ++j) {
	for(int i = minX; i < maxX; ++i) {
		float alpha = ((x[1]*y[2] - x[2]*y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)/ area;
		float beta =  ((x[2]*y[0] - x[0]*y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j) / area;
		float gamma = ((x[0]*y[1] - x[1]*y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j) / area;

		if (alpha >= 0 && beta >= 0 && gamma >= 0) {
			float a = alpha;
			float b = beta;
			float g = gamma;
			float intz = alpha * z[0] + beta * z[1] + gamma * z[2];

			if(intz < state.image_depth[i + j * state.image_width]) {
				state.image_depth[i + j * state.image_width] = intz;

				for(int k = 0; k < state.floats_per_vertex; ++k) {
					float n = 0;
					switch(state.interp_rules[k]) {
						case interp_type::flat: 
							data[k] = v0.data[k];
							break;
						case interp_type::smooth: 
							n = (a / v0.gl_Position[3]) + (b / v1.gl_Position[3]) + (g / v2.gl_Position[3]);

							alpha = a / (n * v[0].gl_Position[3]);
							beta = b / (n * v[1].gl_Position[3]);
							gamma = g / (n* v[2].gl_Position[3]); 
						case interp_type::noperspective: 
							data[k] = (alpha * v0.data[k]) + (beta * v1.data[k]) + (gamma * v2.data[k]);
							break;			 
						default:
							break;
					}
				}
				fragData.data = data;
				state.fragment_shader(fragData, output, state.uniform_data);	
				state.image_color[i + j * state.image_width] = make_pixel( output.output_color[0] * 255, output.output_color[1] * 255, output.output_color[2] * 255);
			}
		}
	}
    }
    delete []v;
}


