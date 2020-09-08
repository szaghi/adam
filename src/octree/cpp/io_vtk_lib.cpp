#include "io_vtk_lib.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <iostream>

namespace pt = boost::property_tree;

void write_vtr(std::vector<float*> fs, int nx, int ny, int nz, float *xg, float *yg, float *zg, std::string filename, std::vector<std::string> var_names) {

    pt::ptree tree, time, cycle, xt, yt, zt;
    std::string format="ascii"; // "binary";

    //time.put_child("DataArray", {});
    //time.put("DataArray", base64_encode((unsigned char const*) "0.00", 4));
    //time.put("DataArray.<xmlattr>.type", "Float64");
    //time.put("DataArray.<xmlattr>.NumberOfTuples", "1");
    //time.put("DataArray.<xmlattr>.Name", "TIME");
    //time.put("DataArray.<xmlattr>.format", "binary");
    //tree.add_child("VTKFile.RectilinearGrid.FieldData.DataArray", time.get_child("DataArray"));

    //cycle.put_child("DataArray", {});
    //cycle.put("DataArray", base64_encode((unsigned char const*) "0001", 4));
    //cycle.put("DataArray.<xmlattr>.type", "Int64");
    //cycle.put("DataArray.<xmlattr>.NumberOfTuples", "1");
    //cycle.put("DataArray.<xmlattr>.Name", "CYCLE");
    //cycle.put("DataArray.<xmlattr>.format", "binary");
    //tree.add_child("VTKFile.RectilinearGrid.FieldData.DataArray", cycle.get_child("DataArray"));

    size_t pad, l;
    unsigned char * tt;

    l = 4*(nx+1); pad = (3-l%3)%3;
    tt = (unsigned char*) malloc(l+pad); memcpy(tt, (unsigned char *) xg, l);
    xt.put_child("DataArray", {});
    std::string ts;
    for(int i=0;i<nx+1;i++) {
        ts += std::to_string(xg[i]);
        if(i<nx) ts += " ";
    }
    xt.put("DataArray", ts.c_str());
    //xt.put("DataArray", b64_encode(tt, l+pad));
    //xt.put("DataArray", base64_encode(tt, l+pad));
    //xt.put("DataArray", base64_encode((unsigned char const*) xg, 4*(nx+1)));
    xt.put("DataArray.<xmlattr>.type", "Float32");
    xt.put("DataArray.<xmlattr>.NumberOfComponents", "1");
    xt.put("DataArray.<xmlattr>.Name", "X");
    xt.put("DataArray.<xmlattr>.format", format);
    tree.add_child("VTKFile.RectilinearGrid.Piece.Coordinates.DataArray", xt.get_child("DataArray"));
    free(tt);

    ts = "";
    for(int i=0;i<ny+1;i++) {
        ts += std::to_string(yg[i]);
        if(i<ny) ts += " ";
    }
    l = 4*(ny+1); pad = (3-l%3)%3;
    tt = (unsigned char*) malloc(l+pad); memcpy(tt, (unsigned char *) yg, l);
    yt.put_child("DataArray", {});
    //yt.put("DataArray", b64_encode(tt, l+pad));
    yt.put("DataArray", ts.c_str());
    yt.put("DataArray.<xmlattr>.type", "Float32");
    yt.put("DataArray.<xmlattr>.NumberOfComponents", "1");
    yt.put("DataArray.<xmlattr>.Name", "Y");
    yt.put("DataArray.<xmlattr>.format", format);
    tree.add_child("VTKFile.RectilinearGrid.Piece.Coordinates.DataArray", yt.get_child("DataArray"));
    free(tt);

    ts = "";
    for(int i=0;i<nz+1;i++) {
        ts += std::to_string(zg[i]);
        if(i<nz) ts += " ";
    }
    l = 4*(nz+1); pad = (3-l%3)%3;
    tt = (unsigned char*) malloc(l+pad); memcpy(tt, (unsigned char *) zg, l);
    zt.put_child("DataArray", {});
    //zt.put("DataArray", b64_encode(tt, l+pad));
    zt.put("DataArray", ts.c_str());
    zt.put("DataArray.<xmlattr>.type", "Float32");
    zt.put("DataArray.<xmlattr>.NumberOfComponents", "1");
    zt.put("DataArray.<xmlattr>.Name", "Z");
    zt.put("DataArray.<xmlattr>.format", format);
    tree.add_child("VTKFile.RectilinearGrid.Piece.Coordinates.DataArray", zt.get_child("DataArray"));
    free(tt);

    std::string extent = "+0 +"+std::to_string(nx)+" +0 +"+std::to_string(ny)+" +0 +"+std::to_string(nz);

    tree.put("VTKFile.RectilinearGrid.<xmlattr>.WholeExtent", extent);
    tree.put("VTKFile.RectilinearGrid.Piece.<xmlattr>.Extent", extent);
    tree.put("VTKFile.<xmlattr>.type","RectilinearGrid");
    tree.put("VTKFile.<xmlattr>.version","1.0");
    tree.put("VTKFile.<xmlattr>.byte_order","LittleEndian");

    int i_var = 0;
    for(auto &f: fs) {
        ts = "";
        for(int k=0;k<nz;k++) 
        for(int j=0;j<ny;j++)
        for(int i=0;i<nx;i++) {
            ts += std::to_string(f[k*nx*ny+j*nx+i]);
            ts += " ";
        }
        pt::ptree temp;
        temp.put_child("DataArray", {});
        temp.put("DataArray", ts.c_str());
        //temp.put("DataArray", base64_encode((unsigned char const*) f, 4*nx*ny*nz));
        temp.put("DataArray.<xmlattr>.type", "Float32");
        temp.put("DataArray.<xmlattr>.NumberOfComponents", "1");
        temp.put("DataArray.<xmlattr>.Name", var_names[i_var]);
        temp.put("DataArray.<xmlattr>.format", format);
        tree.add_child("VTKFile.RectilinearGrid.Piece.CellData.DataArray", temp.get_child("DataArray"));
        i_var++;
    }

    pt::write_xml(filename, tree);

}

void write_vtm(int n_blocks, std::string & filename_all, std::vector<std::string> &filenames) {
    pt::ptree tree;
    int i_block = 0;
    for( auto &fn: filenames) {
        pt::ptree temp;
        temp.put("Block.<xmlattr>.index",i_block);
        std::string bs = std::to_string(i_block+1);
        temp.put("Block.<xmlattr>.name",std::string(5 - bs.length(), '0') + bs);
        //temp.put_child("Block.Dataset", {});
        temp.put("Block.DataSet.<xmlattr>.index", 0);
        temp.put("Block.DataSet.<xmlattr>.file", filenames[i_block]);
        temp.put("Block.DataSet.<xmlattr>.name", filenames[i_block]);
        tree.add_child("VTKFile.vtkMultiBlockDataSet.Block", temp.get_child("Block"));
        i_block++;
    }
    tree.put("VTKFile.<xmlattr>.type", "vtkMultiBlockDataSet"); 
    tree.put("VTKFile.<xmlattr>.version", "1.0");
    tree.put("VTKFile.<xmlattr>.byte_order","LittleEndian");
    pt::write_xml(filename_all, tree);
}

static const std::string base64_chars = 
"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
"abcdefghijklmnopqrstuvwxyz"
"0123456789+/";


static inline bool is_base64(unsigned char c) {
    return (isalnum(c) || (c == '+') || (c == '/'));
}

std::string base64_encode(unsigned char const* bytes_to_encode, unsigned int in_len) {
    std::string ret;
    int i = 0;
    int j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
        char_array_3[i++] = *(bytes_to_encode++);
        if (i == 3) {
            char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
            char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
            char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
            char_array_4[3] = char_array_3[2] & 0x3f;

            for(i = 0; (i <4) ; i++)
                ret += base64_chars[char_array_4[i]];
            i = 0;
        }
    }

    if (i)
    {
        for(j = i; j < 3; j++)
            char_array_3[j] = '\0';

        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for (j = 0; (j < i + 1); j++)
            ret += base64_chars[char_array_4[j]];

        while((i++ < 3))
            ret += '=';

    }

    return ret;

}

std::string base64_decode(std::string const& encoded_string) {
    size_t in_len = encoded_string.size();
    size_t i = 0;
    size_t j = 0;
    int in_ = 0;
    unsigned char char_array_4[4], char_array_3[3];
    std::string ret;

    while (in_len-- && ( encoded_string[in_] != '=') && is_base64(encoded_string[in_])) {
        char_array_4[i++] = encoded_string[in_]; in_++;
        if (i ==4) {
            for (i = 0; i <4; i++)
                char_array_4[i] = static_cast<unsigned char>(base64_chars.find(char_array_4[i]));

            char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
            char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
            char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

            for (i = 0; (i < 3); i++)
                ret += char_array_3[i];
            i = 0;
        }
    }

    if (i) {
        for (j = i; j <4; j++)
            char_array_4[j] = 0;

        for (j = 0; j <4; j++)
            char_array_4[j] = static_cast<unsigned char>(base64_chars.find(char_array_4[j]));

        char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
    }

    return ret;
}

/**
 * `encode.c' - b64
 *
 * copyright (c) 2014 joseph werle
 */

#include <stdio.h>
#include <stdlib.h>
#include "b64.h"

#ifdef b64_USE_CUSTOM_MALLOC
extern void* b64_malloc(size_t);
#endif

#ifdef b64_USE_CUSTOM_REALLOC
extern void* b64_realloc(void*, size_t);
#endif

char *
b64_encode (const unsigned char *src, size_t len) {
    int i = 0;
    int j = 0;
    char *enc = NULL;
    size_t size = 0;
    unsigned char buf[4];
    unsigned char tmp[3];

    // alloc
    enc = (char *) b64_malloc(1);
    if (NULL == enc) { return NULL; }

    // parse until end of source
    while (len--) {
        // read up to 3 bytes at a time into `tmp'
        tmp[i++] = *(src++);

        // if 3 bytes read then encode into `buf'
        if (3 == i) {
            buf[0] = (tmp[0] & 0xfc) >> 2;
            buf[1] = ((tmp[0] & 0x03) << 4) + ((tmp[1] & 0xf0) >> 4);
            buf[2] = ((tmp[1] & 0x0f) << 2) + ((tmp[2] & 0xc0) >> 6);
            buf[3] = tmp[2] & 0x3f;

            // allocate 4 new byts for `enc` and
            // then translate each encoded buffer
            // part by index from the base 64 index table
            // into `enc' unsigned char array
            enc = (char *) b64_realloc(enc, size + 4);
            for (i = 0; i < 4; ++i) {
                enc[size++] = b64_table[buf[i]];
            }

            // reset index
            i = 0;
        }
    }

    // remainder
    if (i > 0) {
        // fill `tmp' with `\0' at most 3 times
        for (j = i; j < 3; ++j) {
            tmp[j] = '\0';
        }

        // perform same codec as above
        buf[0] = (tmp[0] & 0xfc) >> 2;
        buf[1] = ((tmp[0] & 0x03) << 4) + ((tmp[1] & 0xf0) >> 4);
        buf[2] = ((tmp[1] & 0x0f) << 2) + ((tmp[2] & 0xc0) >> 6);
        buf[3] = tmp[2] & 0x3f;

        // perform same write to `enc` with new allocation
        for (j = 0; (j < i + 1); ++j) {
            enc = (char *) b64_realloc(enc, size + 1);
            enc[size++] = b64_table[buf[j]];
        }

        // while there is still a remainder
        // append `=' to `enc'
        while ((i++ < 3)) {
            enc = (char *) b64_realloc(enc, size + 1);
            enc[size++] = '=';
        }
    }

    // Make sure we have enough space to add '\0' character at end.
    enc = (char *) b64_realloc(enc, size + 1);
    enc[size] = '\0';

    return enc;
}


