document.addEventListener('DOMContentLoaded', () => {
    const inputFormula = document.getElementById('input-formula');
    const predictBtn = document.getElementById('predict-btn');

    // 监听输入框点击事件，自动全选文本
    inputFormula.addEventListener('click', () => {
        inputFormula.select();
    });

    // 监听输入框变化事件
    inputFormula.addEventListener('input', () => {
        if (inputFormula.value.trim() !== "") {
            predictBtn.classList.add('active');
        } else {
            predictBtn.classList.remove('active');
        }
    });

    // 监听输入框的 Enter 键事件
    inputFormula.addEventListener('keydown', (event) => {
        if (event.key === 'Enter') {
            event.preventDefault(); // 阻止默认的表单提交行为
            inputFormula.readOnly = true; // 锁定输入框内容
            predictBtn.classList.add('active');
        } else if (event.key === 'Backspace' || event.key === 'Delete') {
            inputFormula.value = ""; // 清空输入框内容
            predictBtn.classList.remove('active');
        }
    });

    predictBtn.addEventListener('click', () => {
        const formula = document.getElementById('input-formula').value;
        const files = document.getElementById('file-upload').files;
        const outputFormat = document.getElementById('output-format').value;
        const outputFormulaChecked = document.getElementById('output-formula').checked;

        console.log('Formula:', formula);
        console.log('Files:', files);
        console.log('Output Format (Digits):', outputFormat);
        console.log('Output Formula:', outputFormulaChecked);

        // 这里可以添加进一步处理逻辑，例如发送请求到服务器进行计算
    });

    function validationSmiles(smiles) {
        // 这里可以添加进一步的 SMILES 格式验证逻辑
        fetch('/validation_smiles', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ smiles: smiles })
        })
        .then(response => response.json())
        .then(data => {
            if (data.valid) {
                inputFormula.readOnly = true;
                predictBtn.classList.add('active');
            } else {
                alert('Invalid SMILES format.');
            }
        })
        .catch(error => {
            console.error('Error:', error);
        });
    }
});

