//============================================================================
// Name        : ahencode.cpp
// Author      : Felix
// Version     : 1.0
// Copyright   : 2016 All rights reserved
// Description : Adaptive Huffman Vitter Invariant Encoder
//============================================================================

#include <iostream>
#include <sstream>
#include <fstream> //for input and output
#include <string>
#include <locale>
#include <cstring>
#include <cmath>
#include "HuffmanTree.h"

using namespace std;

//returns whether node is a left child/not a right child
bool Node::isLeftChild(){
	if(isRoot()){
		return false;
	}

	if(parent->leftChild == this){
		return true;
	} else {
		return false;
	}
}

//initialise the huffman tree with a root node which is the NYT
HuffmanTree::HuffmanTree(){
	rootNode = new Node();
	NYTNode = rootNode;

	//initialise all leaf entries to null (we haven't added any encodings yet!)
	for (int i = 0; i < 256; i++){
		messageToLeaf[i] = NULL;
	}

}


//given a 8-bit number, return the equivalent binary string
std::string getBinary (int message){
	string binaryString = "";

	for (int i = 7; i >= 0; i--){
		if (message >= pow(2, i)) {
			message -= pow(2, i);

			//add encoding to message
			binaryString += "1";
		} else {
			binaryString += "0";
		}

	}

	return binaryString;
}

//given a binary string convert to equivalent decimal value
int getASCII(std::string code){

	int numVal = 0;

	for (int i = code.length()-1; i >= 0; i--){
		if (code[i] == '1') {
			numVal += pow(2, (7-i));
		}

	}

	return numVal;
}



//auxillary function for getting the end of the list
Node* HuffmanTree::getLastListNode(){
	Node* lastNode = rootNode;
	while (lastNode->nextInList != NULL){
		lastNode = lastNode->nextInList;
	}
	return lastNode;
}


//swap with the head of the block of nodes of same weight and type
void HuffmanTree::swapWithHeadBlock(Node* swapNode){

	//if we are the root node we don't do nothing (lol)
	if(swapNode->isRoot()){
		return;
	}

	//get block (nodes of the same weight or weight + 1 if its an internal node)
	int blockWeight = swapNode->weight;

	Node* currNode = rootNode;

	//find node before the block
	while(currNode->nextInList != NULL){
		//case where we have reached our node already, in that case we don't do anything, this is because we do not shift our node down, only up
		if(currNode == swapNode) break;

		//if we are leaf
		if(swapNode->isLeaf() && currNode->nextInList->weight <= blockWeight && currNode->nextInList->isLeaf()) break;

		//if we are internal node
		if(!swapNode->isLeaf() && currNode->nextInList->weight <= blockWeight && !currNode->nextInList->isLeaf()) break;

		currNode = currNode->nextInList;
	}

	//case where it wants to swap with itself
	if(currNode->nextInList == swapNode){
		return;
	}

	Node* swapNodePrev = rootNode;
	//find swapNode's previous
	while(swapNodePrev->nextInList != NULL){

		if(swapNodePrev->nextInList == swapNode) break;

		swapNodePrev = swapNodePrev->nextInList;
	}


	//swap time

	//swap parents
	Node* swapNodeParent = swapNode->parent;
	bool swapNodeIsLeftChild = swapNode->isLeftChild();
	swapNode->parent = currNode->nextInList->parent;
	if(currNode->nextInList->isLeftChild()){
		currNode->nextInList->parent->leftChild = swapNode;
	} else {
		currNode->nextInList->parent->rightChild = swapNode;
	}
	currNode->nextInList->parent = swapNodeParent;
	if(swapNodeIsLeftChild){
		swapNodeParent->leftChild = currNode->nextInList;
	} else {
		swapNodeParent->rightChild = currNode->nextInList;
	}

	//swap in list
	if(currNode->nextInList == swapNodePrev){
		//special case swap with adjacent node in list
		currNode->nextInList = swapNode;
		swapNodePrev->nextInList = swapNode->nextInList;
		swapNode->nextInList = swapNodePrev;
	} else {
		//normal swap
		Node* swapNodeNext = swapNode->nextInList;
		swapNode->nextInList = currNode->nextInList->nextInList;

		swapNodePrev->nextInList = currNode->nextInList;
		currNode->nextInList->nextInList = swapNodeNext;
		currNode->nextInList = swapNode;
	}



}

//add to tree, return the appropriate encoding (for NYT or non-NYT)
std::string HuffmanTree::addToTree(int messageCode){
	string encodedMessage = "";
	
	//NYT Case: if we don't have an entry yet
	if(!alreadyExists(messageCode)){
		//create a two new nodes, make current NYT the parent of both and set left child to be new NYT, right child the new entry
		Node* newLeft = new Node();
		Node* newRight = new Node();
		newLeft->setParent(NYTNode);
		newRight->setParent(NYTNode);
		NYTNode->setLeftChild(newLeft);
		NYTNode->setRightChild(newRight);
		newRight->setMessageCode(messageCode);

		//add to our list of leaves and internal nodes to the end (NYT is internal, newleft after it as is a leaf)
		Node* lastListNode = getLastListNode();
		lastListNode->nextInList = NYTNode;
		NYTNode->nextInList = newRight;

		//update look up array
		if(!setMessageToLeafArray(messageCode, newRight)){
			cout << "Error: couldn't add to array" << endl;
			cout.flush();
		}

		//return the NYT encoding + ASCII
		encodedMessage += getCodeMessage(NYTNode);
		encodedMessage += getBinary(messageCode);

		//update new NYT
		NYTNode = newLeft;

	} else {
		//do nothing
		//we find the leaf node and return that
		encodedMessage = getCodeMessage(getLeafFromMessage(messageCode));
		//note: increment occurs during update tree
	}


	return encodedMessage;
}

//slide and increment
//given a node, performs vitter's slide and increment function
//returns next node to slide and increment
//difference between internal and leaf node handled
Node* HuffmanTree::slideAndIncrement(Node* nodeToInc){

	if (nodeToInc == NULL) return NULL;

	Node* currNode = rootNode;

	//find block to slide and increment (just before the start of block)
	while(currNode->nextInList != NULL){
		//case where we have reached our node already, in that case we don't do anything, this is because we do not shift our node down, only up
		if(currNode == nodeToInc) break;

		//if we are leaf, find start of block of internal nodes (is also the start of block of same weight)
		if(nodeToInc->isLeaf() && currNode->nextInList->weight <= nodeToInc->weight) break; //the point where the block head is finally the same size is less

		//if we are internal node, find start of block of leaf nodes of same weight
		if(!nodeToInc->isLeaf() && currNode->nextInList->weight <= (nodeToInc->weight + 1) && currNode->nextInList->isLeaf()) break;

		currNode = currNode->nextInList;
	}

	//if nodeToInc is already at the head of the block, or is the preblock node area already
	//works for case of root node
	if(currNode->nextInList == nodeToInc || currNode == nodeToInc){

		nodeToInc->setWeight(nodeToInc->weight + 1);

		return nodeToInc->parent;
	}

	//SPECIAL CASE: SWAP CHILD AND PARENT AND PARENT WITHC CHILD
	//This might not happen because the block doesn't include nodes of different type when we swap, once we swap, if we are a leaf we will
	//if we are only swapping with our parent (check that the node after the parent is the child)
	//keep in case
	if(currNode->nextInList == nodeToInc->parent){
		if(nodeToInc->parent->nextInList == nodeToInc){
			//we cannot swap parent with child, and child with parent

			nodeToInc->setWeight(nodeToInc->weight + 1);

			return nodeToInc->parent;
		}
	}

	//couldn't find the block to slide and increment... something went really wrong
	if(currNode->nextInList == NULL) return NULL;


	//process start sliding the entire block downwards
	//save some details about nodeToInc
	Node* nodeToIncPrevParent = nodeToInc->parent;
	bool nodeToIncIsLeftChild = nodeToInc->isLeftChild();
	Node* afterIncNode = nodeToInc->nextInList; //we will break this off our list so we have to keep track to connect to IncNode's previous
	//save the previous start of the block
	Node* nextNode = currNode->nextInList;

	//place our node at the start of the block
	currNode->nextInList = nodeToInc;
	nodeToInc->nextInList = nextNode;

	//now that we've placed our node in the right position, we can start going through and updating the nodes to their right position in the block
	currNode = currNode->nextInList;

	//while we are still in the block and not at the list's end
	while(currNode->nextInList != NULL){
		if(nodeToInc->isLeaf() && currNode->weight != nodeToInc->weight) break;

		//SPECIAL CASE
		//for internal node we may reach ourselves, but we are not a member of the block yet (the block of weight+1)
		//so we need to deal with case where we are an internal node, not in the block but we are also the node to swap...
		//if we are the node to swap then we must leave loop
		if(!nodeToInc->isLeaf() && currNode->weight != (nodeToInc->weight+1)) {
			if(currNode != nodeToInc){
				break;
			}
		}


		//case where we've reached our node's original position
		//if we have reached nodeToInc's previous, now we want to connect to the node in front of nodeToInc
		//and give it nodeToInc's previous parents
		if (currNode->nextInList == nodeToInc){
			currNode->nextInList = afterIncNode;

			//update this parent to nodeToInc's parent
			if(nodeToIncIsLeftChild){
				nodeToIncPrevParent->leftChild = currNode;
			} else {
				nodeToIncPrevParent->rightChild = currNode;
			}
			currNode->parent = nodeToIncPrevParent;

			//now we stop because we've finished shifting!
			break;
		} else {

			//update this parents with the node in front of our node
			if(nextNode->isLeftChild()){
				nextNode->parent->leftChild = currNode;
			} else {
				nextNode->parent->rightChild = currNode;
			}
			currNode->parent = nextNode->parent;
		}

		//iterate through
		currNode = currNode->nextInList;
		nextNode = currNode->nextInList;
	}

	//increment node
	nodeToInc->setWeight(nodeToInc->weight + 1);

	//return the next node to slide and increment
	Node* nextNodeToSlideAndIncrement = NULL;
	if(nodeToInc->isLeaf()){
		//leaf nodes return its new parent
		nextNodeToSlideAndIncrement = nodeToInc->parent;
	} else {
		//internal nodes return its old parent (before slide and increment)
		nextNodeToSlideAndIncrement = nodeToIncPrevParent;
	}

	return nextNodeToSlideAndIncrement;
}


//update stage after adding message to tree
void HuffmanTree::updateTree(Node* currNode){
	if(currNode == NULL) return;

	//we swap our node with the head of the block
	swapWithHeadBlock(currNode);

	//slide and increment
	currNode = slideAndIncrement(currNode);
	while (currNode != NULL){
		currNode = slideAndIncrement(currNode);
	}
}


//given a code decode into original message
std::string HuffmanTree::getCodeMessage(Node* leafNode){
	if(leafNode == NULL || this == NULL) return "";
	//iterate up tree, appending to message
	std::string reversedCode = ""; //reversed as we are going backwards up the tree
	Node* currNode = leafNode;

	while(!currNode->isRoot()){

		if(currNode->isLeftChild()){
			reversedCode += HuffmanTree::LEFT_CHILD_BINARY;
		} else {
			reversedCode += HuffmanTree::RIGHT_CHILD_BINARY;
		}

		currNode = currNode->parent;
	}

	std::string codeMessage = "";

	//reverse the code
	for(int n = reversedCode.length()-1; n >= 0; n--){
		codeMessage += reversedCode[n];
	}

	return codeMessage;
}

//encode the given message and return the adaptive huffman code
std::string HuffmanTree::encodeMessageCode(int messageCode){
	std::string encodedMessage = addToTree(messageCode);

	//update tree (swap with head of block, then slide and increment through block)
	updateTree(getLeafFromMessage(messageCode));

	return encodedMessage;
}

std::string HuffmanTree::encodeMessage(std::wstring inputLine, bool spaceDelimSet){
	std::string completeEncodedMessage  = "";
	for(unsigned int i = 0; i < inputLine.length(); i++){	
		//(char)inputLine[i] is different to (int) inputLine[i]
		//difference between 227 and -29
		//because (char) is an signed int... holy crap
		//char treats it like an ASCII value, so anything above 128 goes negative
		int numVal = (int) inputLine[i];
		string encodedMessage = encodeMessageCode(numVal);
		completeEncodedMessage += encodedMessage;

		if(spaceDelimSet && i < (inputLine.length() - 1)){
			completeEncodedMessage += " ";
		}

	}

	return completeEncodedMessage;
}


int main(int argc, char* argv[]){
	//set to UTF-8 for linux
	setlocale(LC_ALL, "");
	//read any arguments
	bool spaceDeliminterSet = false;
	if(argc > 1){
		string argument1 = argv[1];

		if(argument1.compare("-s") == 0){
			spaceDeliminterSet = true;
		}
	}

	//input
	for (wstring inputLine; getline(std::wcin, inputLine);) {
		
		HuffmanTree* huffmanTree = new HuffmanTree();
		std::string encodedMessage = huffmanTree->encodeMessage(inputLine, spaceDeliminterSet);
		cout << encodedMessage  << endl;

	}

	return 0;
}

