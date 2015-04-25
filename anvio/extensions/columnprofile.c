#include <Python.h>
#include <stdio.h>
#include <string.h>
#include "structmember.h"


typedef struct {
    PyObject_HEAD
    PyObject *column;
    int coverage;
    int pos;
    PyObject *split_name;
    int sample_id;
    PyObject *test_class;
    PyObject *profile;
} ColumnProfile;

static void
ColumnProfile_dealloc(ColumnProfile* self)
{

    Py_XDECREF(self->column);
    Py_XDECREF(self->split_name);
    Py_XDECREF(self->test_class);
    Py_XDECREF(self->profile);

    self->ob_type->tp_free((PyObject*)self);
}


static PyObject *
ColumnProfile_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    ColumnProfile *self;

    self = (ColumnProfile *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->column = PyString_FromString("");
        if (self->column == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        
        self->split_name = PyString_FromString("");
        if (self->split_name == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }
        
        self->profile = PyDict_New();
        if (self->profile == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }

        PyObject *variability_module = PyImport_ImportModule("anvio.variability");
        PyObject *variablity_class = PyObject_GetAttrString(variability_module, "VariablityTestFactory");
        self->test_class = PyInstance_NewRaw(variablity_class, NULL);     

        if (self->test_class == NULL)
          {
            Py_DECREF(self);
            return NULL;
          }

        self->coverage = 0;
        self->pos = 0;
        self->sample_id = 0;
    }

    return (PyObject *)self;
}


int count_chars(const char* string, char ch)
{
    int count = 0;
    for(; *string; count += (*string++ == ch));
    return count;
}

int compare ( const void *pa, const void *pb )
{
    const int *a = pa;
    const int *b = pb;

    if(a[0] == b[0])
        return a[1] - b[1];
    else
        return a[0] - b[0];
}


static int
ColumnProfile_init(ColumnProfile *self, PyObject *args, PyObject *kwds)
{
    PyObject *column=NULL, *split_name=NULL, *test_class=Py_None, *tmp;

    static char *kwlist[] = {"column", "coverage", "pos", "split_name", "sample_id", "test_class", NULL};

    if (! PyArg_ParseTupleAndKeywords(args, kwds, "|OiiOiO", kwlist, 
                                      &column,
                                      &self->coverage,
                                      &self->pos,
                                      &split_name, 
                                      &self->sample_id,
                                      &test_class))
        return -1; 

    if (column) {
        tmp = self->column;
        Py_INCREF(column);
        self->column = column;
        Py_XDECREF(tmp);
    }

    if (split_name) {
        tmp = self->split_name;
        Py_INCREF(split_name);
        self->split_name = split_name;
        Py_XDECREF(tmp);
    }

    if (test_class) {
        tmp = self->test_class;
        Py_INCREF(test_class);
        self->test_class = test_class;
        Py_XDECREF(tmp);
    }

    PyDict_SetItemString(self->profile, "sample_id", PyInt_FromLong(self->sample_id));
    PyDict_SetItemString(self->profile, "split_name", self->split_name);
    PyDict_SetItemString(self->profile, "pos", PyInt_FromLong(self->pos));
    PyDict_SetItemString(self->profile, "consensus", Py_None);
    PyDict_SetItemString(self->profile, "coverage", (self->coverage != 0) ? PyInt_FromLong(self->coverage) : PyInt_FromLong(PyString_Size(self->column)));
    PyDict_SetItemString(self->profile, "n2n1ratio", PyFloat_FromDouble(0));
    PyDict_SetItemString(self->profile, "competing_nts", Py_None);

    const char nucleotides[6] = "ATCGN";
    const char * _column = PyString_AsString(self->column);

    int i;
    int nucleotides_arr[5][2];
    char np[2] = "";
    char competing_nts[3] = "";
    double n2n1ratio;
    double minratio;

    for (i=0; i<5; i++){
        int cc = count_chars(_column, nucleotides[i]); 
        nucleotides_arr[i][0] = cc;    
        nucleotides_arr[i][1] = i;
        np[0] = nucleotides[(int)nucleotides_arr[i][1]];
        PyDict_SetItemString(self->profile, np, PyInt_FromLong(cc));
    }

    qsort(nucleotides_arr, 5, sizeof nucleotides_arr[0], compare);

    np[0] = nucleotides[(int)nucleotides_arr[4][1]];

    if (!nucleotides_arr[3][0]) {
        PyDict_SetItemString(self->profile, "consensus", PyString_FromString(np));
        return 0;
    }

    competing_nts[0] = nucleotides[(int)nucleotides_arr[4][1]];
    competing_nts[1] = nucleotides[(int)nucleotides_arr[3][1]];

    if (nucleotides_arr[4][1] == 4) {
        return 0;
    }

    PyDict_SetItemString(self->profile, "consensus", PyString_FromString(np));

    if (nucleotides_arr[3][1] == 4) {
        return 0;
    }

    n2n1ratio = ((double)nucleotides_arr[3][0] / (double)nucleotides_arr[4][0]);
    PyObject *min_acceptable_ratio_given_coverage = PyObject_CallMethod(self->test_class, "min_acceptable_ratio_given_coverage", "I", self->coverage);
    minratio = PyFloat_AsDouble(min_acceptable_ratio_given_coverage);

    if (test_class != Py_None) {
        if (n2n1ratio > minratio){
            PyDict_SetItemString(self->profile, "n2n1ratio", PyFloat_FromDouble(n2n1ratio));
            PyDict_SetItemString(self->profile, "competing_nts", PyString_FromString(competing_nts));        
        }
    }

    else {
        PyDict_SetItemString(self->profile, "n2n1ratio", PyFloat_FromDouble(n2n1ratio));
        PyDict_SetItemString(self->profile, "competing_nts", PyString_FromString(competing_nts));   
    }

    return 0;
}


static PyMemberDef ColumnProfile_members[] = {
    {"column", T_OBJECT_EX, offsetof(ColumnProfile, column), 0,
     "column"},
    {"coverage", T_INT, offsetof(ColumnProfile, coverage), 0,
     "coverage"},
    {"pos", T_INT, offsetof(ColumnProfile, pos), 0,
     "pos"},          
    {"split_name", T_OBJECT_EX, offsetof(ColumnProfile, split_name), 0,
     "split_name"},
    {"sample_id", T_INT, offsetof(ColumnProfile, sample_id), 0,
     "sample_id"},
    {"test_class", T_OBJECT_EX, offsetof(ColumnProfile, test_class), 0,
     "test_class"},          
    {"profile", T_OBJECT_EX, offsetof(ColumnProfile, profile), 0,
     "profile"},        
    {NULL}  /* Sentinel */
};

static PyTypeObject ColumnProfileType = {
    PyObject_HEAD_INIT(NULL)
    0,                                          /* ob_size */
    "columnprofile.ColumnProfile",              /* tp_name */
    sizeof(ColumnProfile),                      /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)ColumnProfile_dealloc,          /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_compare */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "ColumnProfile objects",                    /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    0,                                          /* tp_methods */
    ColumnProfile_members,                      /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)ColumnProfile_init,               /* tp_init */
    0,                                          /* tp_alloc */
    ColumnProfile_new,                          /* tp_new */
};

static PyMethodDef module_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initcolumnprofile(void) 
{
    PyObject* m;

    if (PyType_Ready(&ColumnProfileType) < 0)
        return;

    m = Py_InitModule3("columnprofile", module_methods,
                       "Example ColumnProfile class.");

    if (m == NULL)
      return;

    Py_INCREF(&ColumnProfileType);
    PyModule_AddObject(m, "ColumnProfile", (PyObject *)&ColumnProfileType);
}
