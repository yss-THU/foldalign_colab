#ifndef STACK_SSL
#define STACK_SSL

#include <iostream>

/******************************************************************************
*                                                                             *
*   Copyright 2004 - 2007 Jakob Hull Havgaard, hull@bioinf.ku.dk              *
*                                                                             *
*   This file is part of Foldalign                                            *
*                                                                             *
*   Foldalign is free software; you can redistribute it and/or modify         *
*   it under the terms of the GNU General Public License as published by      *
*   the Free Software Foundation; either version 2 of the License, or         *
*   (at your option) any later version.                                       *
*                                                                             *
*   Foldalign is distributed in the hope that it will be useful,              *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*   GNU General Public License for more details.                              *
*                                                                             *
*   You should have received a copy of the GNU General Public License         *
*   along with Foldalign; if not, write to the Free Software                  *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA *
*                                                                             *
******************************************************************************/

template< class item, class index >
class stack_ssl {

public:
	stack_ssl() : current(0), num(0) {pthread_mutex_init(&mutex, NULL);};
	
	void push(const item value) {
		pthread_mutex_lock(&mutex);
		node* new_node = new node();
		new_node->value = value;
		new_node->previous = current;
		current = new_node;
		num++;
		pthread_mutex_unlock(&mutex);
	}
	
	item* push_ptr(const item value) {
		pthread_mutex_lock(&mutex);
		node* new_node = new node();
		new_node->value = value;
		new_node->previous = current;
		current = new_node;
		num++;
		item* returnValue= &(new_node->value);
		pthread_mutex_unlock(&mutex);
		return returnValue;
	}
	
	
	item pop();
	
	index size() const {return num;};
	
	~stack_ssl() {
	
		while (current != 0) {
			node* pre = current->previous;
			delete current;
			current = pre;
		}
		pthread_mutex_destroy(&mutex);

	};
	
private:

	struct node {
		item value;
		node* previous;
	};
	
	node* current;

	index num;
	pthread_mutex_t mutex;
};

template< class item, class index >
inline item stack_ssl< item, index>::pop() {

	pthread_mutex_lock(&mutex);
	if (num == 0) {
		std::cerr << "Program error! Stack underflow" << std::endl;
		pthread_mutex_unlock(&mutex);
		throw;
	}
	
	item return_value = current->value;
	
	node* old_node = current;
	
	current = current->previous;
	delete old_node;
	
	num--;
	pthread_mutex_unlock(&mutex);

	return return_value;
	
}


#endif /* STACK_SSL */
