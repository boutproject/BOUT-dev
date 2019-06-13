#include <boutcomm.hxx>

Request SingleProcessorCommunicator::irecvBytes(void *data, int length, int destination, int identifier) {
  ASSERT1(destination == 0);

  // Get the request ID
  request = next_request_id++;

  // Insert into map
  posted_requests[request] = {identifier, data, length, false};
  
  return request;
}

void SingleProcessorCommunicator::sendBytes(void *data, int length, int destination, int identifier) {
  ASSERT1(destination == 0);
  
}

int SingleProcessorCommunicator::waitAny(const std::vector<Request> &requests) {
  
}

