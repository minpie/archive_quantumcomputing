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

def OperatorGet_Q(op1):
    # op1의 값을 반환하는 함수.
    '''
    * op1 비트 배치:
    LSB ... MSB
    * 반환값 비트 배치:
    MSB ... LSB
    '''
    result = 0
    for i in range(len(op1)-1, -1, -1):
        result = result << 1
        result = result + int(op1[i])
    return result 

def OperatorSet_QI(op1, op2):
    # op1의 값을 op2로 설정하는 함수.
    '''
    * op1 비트 배치:
    LSB ... MSB
    * op2 비트 배치:
    MSB ... LSB
    '''
    for i in range(0, len(op1)):
        bit = (op2 >> i) & 0b1
        if(bit == 1):
            X | op1[i]
    return

def GetBitLen(num1):
    result = 0
    if(num1 < 2):
        result = 1
    else:
        result = int(math.floor(math.log2(num1))) + 1
    return result

def Example1():
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