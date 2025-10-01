# import:
import qiskit
from qiskit import QuantumCircuit
from qiskit.circuit import QuantumRegister, ClassicalRegister
from qiskit_aer import AerSimulator
import math

def OperatorGet_Q(op1):
    # op1의 값을 반환하는 함수.
    '''
    OperatorGet_Q(op1)
    Do:
    - (get value of op1)
    Operands type:
    - op1: QuantumRegister
    Description:
    - this function gets value of op1
    - Work in progress
    '''
    #
    return

def OperatorSet_QI(op1, op2, targetCircuit):
    # op1의 값을 op2로 설정하는 함수.
    '''
    OperatorSet_QI(op1, op2, targetCircuit)
    Do:
    - op1 = op2
    on targetCircuit
    Operands type:
    - op1: QuantumRegister
    - op2: integer
    - targetCircuit: QuantumCircuit
    Description:
    - this function sets op1 to op2 if op1 are zerorized state.
    '''
    #
    for i in range(0, len(op1)):
        bit = (op2 >> i) & 0b1
        if(bit == 1):
            targetCircuit.x(op1[i])
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

    # parameter:
    c_n1 = 234 # integer: n1
    c_sizeofn = 8

    # elements:
    q_n1 = QuantumRegister(c_sizeofn, "n1")
    c_result = ClassicalRegister(c_sizeofn, "result")
    mainCircuit = QuantumCircuit(q_n1, c_result)

    # assignment: n1
    OperatorSet_QI(q_n1, c_n1, mainCircuit)
    #

    # measure:
    mainCircuit.measure(q_n1, c_result)

    # end circuit
    simulator = AerSimulator()
    job = simulator.run(mainCircuit, shots=1)
    result = job.result()


    # analyze result:
    raw = list(result.get_counts(mainCircuit).keys())[0]
    print("c_n1 = {} | q_n1 = {}".format(c_n1, int(raw, 2)))
    pass

def main():
    Example1()
    pass

#### main() entry
if __name__ == "__main__":
    main()
    pass