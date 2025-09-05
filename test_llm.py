import requests
import base64
import os
from openai import OpenAI


"""
deepseek-r1
deepseek-v3
deepseek-v3.1
deepseek-r1-distill-qwen-1.5b
deepseek-r1-distill-qwen-7b
deepseek-r1-distill-qwen-14b
deepseek-r1-distill-qwen-32b
deepseek-r1-distill-llama-8b
deepseek-r1-distill-llama-70b
"""


# 若没有配置环境变量，请用百炼API Key将下行替换为：api_key="sk-xxx",
#  os.getenv("DASHSCOPE_API_KEY"),  # 如何获取API Key：https://help.aliyun.com/zh/model-studio/developer-reference/get-api-key
client = OpenAI(
    api_key="sk-f4207c24fda4439f8d8fd711da49cf8c",
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
)

if __name__ == "__main__":
    # Process an image example
    
    API_KEY = "sk-f4207c24fda4439f8d8fd711da49cf8c"
    
    # 此处以 deepseek-r1 为例，可按需更换模型名称。
    completion = client.chat.completions.create(
        model="deepseek-r1-distill-qwen-7b",
        messages=[{"role": "user", "content": "9.9和9.11谁大"}],
    )
    print("最终答案：")
    print(completion.choices[0].message.content)


