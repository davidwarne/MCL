/* MCL: Monte Carlo Library.
 * Copyright (C) 2017  David J. Warne
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "mcl.h"

/**
 * @brief copies the memory structure of the dataset d
 * @details This function copys all fields and allocates new data arrays based 
 * on field sizes. That is the a new data set structure is create using d as a
 * template.
 *
 * @param d data structure to copy
 * @returns a data structure with same field inforamtion and array sizes as d.
 *
 *
 * @note This function does not populate the new data arrays.
 * @todo add a flag option to copy with data.
 */
Dataset* 
copyDataset(Dataset * d)
{
    size_t i;
    Dataset *s;
    s = (Dataset *)malloc(sizeof(Dataset));
    s->numFields = d->numFields;
    s->fields = (field *)malloc(s->numFields*sizeof(field));
    /*allocate  new data arrays*/
    for (i=0;i<s->numFields;i++)
    {
        s->fields[i] = d->fields[i];
        s->fields[i].data_array = (void *)malloc(d->fields[i].numBytes);
    }
    return s;
}
