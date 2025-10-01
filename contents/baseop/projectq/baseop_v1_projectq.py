# import:
from projectq import MainEngine
from projectq.backends import ClassicalSimulator, ResourceCounter
from projectq.ops import (
    Measure,
    H,
    All,
    Toffoli,
    CNOT,
    X,
    Swap,
)
from projectq.meta import Loop, Compute, Uncompute, Control
import math

# function:
def OperatorGet_Q(op1):
    # op1의 값을 반환하는 함수.
    '''
    OperatorGet_Q(op1)
    Do:
    - (get value of op1)
    Operands type:
    - op1: qubits
    Description:
    - this function gets value of op1
    '''
    #
    result = 0
    for i in range(len(op1)-1, -1, -1):
        result = result << 1
        result = result + int(op1[i])
    return result 

def OperatorSet_QI(op1, op2):
    # op1의 값을 op2로 설정하는 함수.
    '''
    OperatorSet_QI(op1, op2)
    Do:
    - op1 = op2
    Operands type:
    - op1: qubits
    - op2: integer
    Description:
    - this function sets op1 to op2 if op1 are zerorized state.
    '''
    #
    for i in range(0, len(op1)):
        bit = (op2 >> i) & 0b1
        if(bit == 1):
            X | op1[i]
    return

def GetBitLen(op1):
    # op1이 2진수로 몇 비트 크기의 수인지 측정하는 함수.
    '''
    GetBitLen(op1)
    Do:
    - (get length of op1 as bit)
    Operands type:
    - op1: integer
    Description:
    - this function measures the minimum length of op1 as bit.
    - negative value and non-integer value is not considered.
    '''
    #
    result = 0
    if(op1 < 2):
        result = 1
    else:
        result = int(math.floor(math.log2(op1))) + 1
    return result

def Example1():
    '''
    Example circuit 1 shows how use the functions.
    
    '''
    #
    objEngine = MainEngine(ClassicalSimulator())
    
    c_a = 1234
    c_b = 255

    q_a = objEngine.allocate_qureg(GetBitLen(c_a))
    q_b = objEngine.allocate_qureg(GetBitLen(c_b))

    OperatorSet_QI(q_a, c_a)
    OperatorSet_QI(q_b, c_b)

    All(Measure) | q_a
    All(Measure) | q_b

    objEngine.flush()

    print("q_a = {}".format(OperatorGet_Q(q_a)))
    print("q_b = {}".format(OperatorGet_Q(q_b)))
    pass

def main():
    Example1()
    pass

#### main() entry
if __name__ == "__main__":
    main()
    pass